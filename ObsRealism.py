#!/bin/env python

'''
The statistical observational realism suite presented in Bottrell et al 2019b. If you use this suite in your research or for any other purpose, I would appreciate a citation. If you have questions or suggestions on how to improve or broaden the suite, please contact me at cbottrel "at" uvic "dot" ca.

Version 0.3

Update History:
(v0_1) - February-2016 Correction to the way that the poisson noise is handled. Changed incorrect float padding to a sampling of a true Poisson distribution with correct implementation of the Gain quantity. SkyServer 'sky' and 'skysig' quantities are added to the header keywords when using real SDSS images.
(v0_2) - January-2019 Spec2SDSS_gri now incorporates the redshift (factor of (1+z)**5) to the wavelength specific intensities when generating the bandpass AB surface brightnesses. The redshift factor is now removed from intensity scaling step in this version.
(v0_3) - January-2019 Computes ra,dec from the source mask when determining image position. This avoids image registration offsets in row and column numbers for each band. Prepared for public release.
'''

import numpy as np
import os,sys,string,time
from scipy.interpolate import RectBivariateSpline
import scipy.ndimage
import warnings
from astropy.io import fits
from astropy.utils.exceptions import AstropyWarning
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.cosmology import FlatLambdaCDM
import sep
sep.set_extract_pixstack(9999999)

realsim_dir = os.path.dirname(os.path.abspath(__file__))

def rebin(array, dimensions=None, scale=None):
    """
    Return the array 'array' to the new 'dimensions'
    conserving the flux in the bins. The sum of the
    array will remain the same as the original array.
    Congrid from the scipy recipies does not generally
    conserve surface brightness with reasonable accuracy.
    As such, even accounting for the ratio of old and new
    image areas, it does not conserve flux. This function
    nicely solves the problem and more accurately
    redistributes the flux in the output. This function
    conserves FLUX so input arrays in surface brightness
    will need to use the ratio of the input and output
    image areas to go back to surface brightness units.
    
    EXAMPLE
    -------
    
    In [0]:
    
    # input (1,4) array (sum of 6)
    y = np.array([0,2,1,3]).reshape(1,4).astype(float)
    # rebin to (1,3) array
    yy = rebin(y,dimensions=(1,3))
    print yy
    print np.sum(yy)
    
    Out [0]:
    
    Rebinning to Dimensions: 1, 3
    [[0.66666667 2.         3.33333333]]
    6.0

    RAISES
    ------
    AssertionError
        If the totals of the input and result array don't
        agree, raise an error because computation may have
        gone wrong.
        
    Copyright: Martyn Bristow (2015) and licensed under GPL v3:
    i.e. free to use/edit but no warranty.
    """
    if dimensions is not None:
        if isinstance(dimensions, float):
            dimensions = [int(dimensions)] * len(array.shape)
        elif isinstance(dimensions, int):
            dimensions = [dimensions] * len(array.shape)
        elif len(dimensions) != len(array.shape):
            raise RuntimeError('')
    elif scale is not None:
        if isinstance(scale, float) or isinstance(scale, int):
            dimensions = map(int, map(round, map(lambda x: x*scale, array.shape)))
        elif len(scale) != len(array.shape):
            raise RuntimeError('')
    else:
        raise RuntimeError('Incorrect parameters to rebin.\n\trebin(array, dimensions=(x,y))\n\trebin(array, scale=a')
    #print "Rebinning to Dimensions: %s, %s" % tuple(dimensions)
    import itertools
    dY, dX = map(divmod, map(float, array.shape), dimensions)
 
    result = np.zeros(dimensions)
    for j, i in itertools.product(*map(range, array.shape)):
        (J, dj), (I, di) = divmod(j*dimensions[0], array.shape[0]), divmod(i*dimensions[1], array.shape[1])
        (J1, dj1), (I1, di1) = divmod(j+1, array.shape[0]/float(dimensions[0])), divmod(i+1, array.shape[1]/float(dimensions[1]))
        
        # Moving to new bin
        # Is this a discrete bin?
        dx,dy=0,0
        if (I1-I == 0) | ((I1-I == 1) & (di1==0)):
            dx = 1
        else:
            dx=1-di1
        if (J1-J == 0) | ((J1-J == 1) & (dj1==0)):
            dy=1
        else:
            dy=1-dj1
        # Prevent it from allocating outide the array
        I_=min(dimensions[1]-1,I+1)
        J_=min(dimensions[0]-1,J+1)
        result[J, I] += array[j,i]*dx*dy
        result[J_, I] += array[j,i]*(1-dy)*dx
        result[J, I_] += array[j,i]*dy*(1-dx)
        result[J_, I_] += array[j,i]*(1-dx)*(1-dy)
    allowError = 0.1
    assert (abs(array.sum()) < abs(result.sum()) * (1+allowError)) & (abs(array.sum()) > abs(result.sum()) * (1-allowError))
    return result
    
def ObsRealism(inputName,outputName,band='r',
                cosmo=FlatLambdaCDM(H0=70,Om0=0.3),
                common_args = { 
                                'redshift'      : 0.1,   # mock observation redshift
                                'rebin_to_CCD'  : False, # rebin to CCD angular scale
                                'CCD_scale'     : 0.396, # CCD angular scale in [arcsec/pixel]
                                'add_false_sky' : False, # add gaussian sky
                                'false_sky_sig' : 24.2,  # gaussian sky standard dev [AB mag/arcsec2]
                                'add_false_psf' : False, # convolve with gaussian psf
                                'false_psf_fwhm': 1.0,   # gaussian psf FWHM [arcsec]
                                'add_poisson'   : False, # add poisson noise to galaxy
                                'add_sdss_sky'  : False, # insert into real SDSS sky (using sdss_args)
                                'add_sdss_psf'  : False, # convolve with real SDSS psf (using sdss_args)

                              },
               sdss_args    = {
                                'sdss_run'      : 745,       # sdss run
                                'sdss_rerun'    : 40,        # sdss rerun
                                'sdss_camcol'   : 1,         # sdss camcol
                                'sdss_field'    : 517,       # sdss field
                                'sdss_ra'       : 236.1900,  # ra for image centroid
                                'sdss_dec'      : -0.9200,   # ec for image centroid
                              }
               ):
    
    '''
    Add realism to idealized unscaled image.
    
    "redshift": The redshift at which the synthetic image is to be mock-observed. Given that the image should be in surface brightness units and appropriately dimmed by (1+z)^-5, the redshift is only used to determine the angular-to-physical scale of the image -- to which it is appropriately rebinned corresponding to the desired CCD pixel scale.
    
    "rebin_to_CCD": If TRUE, the image is rebinned to the CCD scale identified by the "CCD_scale" keyword. The rebinning is determined by first computing the physical-to-angular scale associated with the target redshift [kpc/arcsec]. Combining this number with the scale of the original image in physical units [kpc/pixel], we obtain the rebinning factor that is neccesary to bring the image to the desired CCD pixel scale [arcsec/pixel].
    
    "CCD_scale": The CCD scale to which the images are rebinned if rebin_to_CCD is TRUE.
    
    "add_false_sky": If TRUE, a Gaussian sky is added to the image with a noise level that is idenfitied by the "false_sky_sig" keyword.
    
    "false_sky_sig": The standard deviation of Gaussian sky that is added to the image if "add_false_sky" is TRUE. The value must be expressed in relative magnitude units (AB mag/arcsec2).
    
    "add_false_psf": If TRUE, a Gaussian PSF is added to the image with a FWHM that is idenfitied by the "false_psf_fwhm" keyword.
    
    "false_psf_fwhm": The FWHM of the PSF that is convolved with the image if "add_false_psf" is TRUE. The value must be expressed in arcsec.
    
    "add_poisson": If TRUE, add Poisson noise to the image using either the calibration info and gain from the real image properties ("add_sdss_sky"=TRUE) or generic values derived from averages over SDSS fields.
    
    "add_sdss_sky": If True, insert into real SDSS sky using arguments in "sdss_args".
    
    "add_sdss_psf": If True and "add_sdss_sky"=True, reconstruct the PSF at the injection location and convolve with the image.
    '''
    
    # mock observation redshift
    redshift = common_args['redshift']
    # speed of light [m/s]
    speed_of_light = 2.99792458e8
    # kiloparsec per arcsecond scale
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z=redshift).value/60. # [kpc/arcsec]
    # luminosity distance in Mpc
    luminosity_distance = cosmo.luminosity_distance(z=redshift) # [Mpc]
    
    # img header and data
    with fits.open(inputName,mode='readonly') as hdul:
        # img header
        header = hdul[0].header
        # img data
        img_data = hdul[0].data
    
#    # header properties
#    sim_tag = header['SIMTAG']
#    sub_tag = header['SUBTAG']
#    isnap = header['ISNAP']
#    axis = header['CAMERA']
#    band = header['FILTER'][0]
#
#    # unique simulID
#    simulID = '{}-{}-{}-{}'.format(sim_tag,sub_tag,isnap,axis)
#
#    band = header['FILTER'][0]

    # collect physical pixel scale
    kpc_per_pixel = header['CDELT1']/1000. # [kpc/pixel]
    # compute angular pixel scale from cosmology
    arcsec_per_pixel = kpc_per_pixel / kpc_per_arcsec # [arcsec/pixel]
     
    # img in AB nanomaggies per arcsec2
    img_nanomaggies = 10**(-0.4*(img_data-22.5)) # [nmgys/arcsec2]
    # apply pixel scale [arcsec/pixel]2 to convert to calibrated flux
    img_nanomaggies *= arcsec_per_pixel**2 # [nmgs]
    # update units of image header to linear calibrated scale
    header['BUNIT'] = 'AB nanomaggies'
    
#    print('\nRaw image:')
#    print('kpc_per_arcsec: {}'.format(kpc_per_arcsec))
#    print('kpc_per_pixel: {}'.format(kpc_per_pixel))
#    print('arcsec_per_pixel: {}'.format(arcsec_per_pixel))
#    m_AB = -2.5*np.log10(np.sum(img_nanomaggies))+22.5
#    print('AB_magnitude: {} at z={}'.format(m_AB,redshift))
#    M_AB = m_AB-5*np.log10(luminosity_distance.value)-25
#    print('AB_Magnitude: {}'.format(M_AB))

    # Add levels of realism
    
    if common_args['rebin_to_CCD']:
        '''
        Rebin image to a given angular CCD scale
        '''
        # telescope ccd angular scale
        ccd_scale = common_args['CCD_scale']
        # axes of original image
        nPixelsOld = img_nanomaggies.shape[0]
        # axes of regridded image
        nPixelsNew = int(np.floor((arcsec_per_pixel/ccd_scale)*nPixelsOld))
        # rebin to new ccd scale
        if nPixelsNew>nPixelsOld:
            interp = RectBivariateSpline(np.linspace(-1,1,nPixelsOld),np.linspace(-1,1,nPixelsOld),img_nanomaggies,kx=1,ky=1)
            img_nanomaggies = interp(np.linspace(-1,1,nPixelsNew),np.linspace(-1,1,nPixelsNew))*(nPixelsOld/nPixelsNew)**2
        else:
            img_nanomaggies = rebin(img_nanomaggies,(nPixelsNew,nPixelsNew))
        # new kpc_per_pixel on ccd
        kpc_per_pixel = kpc_per_arcsec * ccd_scale
        # new arcsec per pixel
        arcsec_per_pixel = ccd_scale
        # header updates
        if nPixelsNew%2: CRPIX = float(nPixelsNew/2)
        else: CRPIX = float(nPixelsNew/2)+0.5
        header['CRPIX1'] = CRPIX
        header['CRPIX2'] = CRPIX
        header['CDELT1'] = kpc_per_pixel*1000
        header['CDELT2'] = kpc_per_pixel*1000
#        print('\nAfter CCD scaling:')
#        print('kpc_per_arcsec: {}'.format(kpc_per_arcsec))
#        print('kpc_per_pixel: {}'.format(kpc_per_pixel))
#        print('arcsec_per_pixel: {}'.format(arcsec_per_pixel))
#        m_AB = -2.5*np.log10(np.sum(img_nanomaggies))+22.5
#        print('AB_magnitude: {} at z={}'.format(m_AB,redshift))
#        M_AB = m_AB-5*np.log10(luminosity_distance.value)-25
#        print('AB_Magnitude: {}'.format(M_AB))

    # convolve with gaussian psf
    if common_args['add_false_psf']:
        '''
        Add Gaussian PSF to image with provided FWHM in
        arcseconds.
        '''
        std = common_args['false_psf_fwhm']/arcsec_per_pixel/2.355
        kernel = Gaussian2DKernel(stddev=std)
        img_nanomaggies = convolve(img_nanomaggies, kernel)
        
    # add poisson noise to image
    if common_args['add_poisson'] and not common_args['add_sdss_sky']:
        '''
        Add shot noise to image assuming the average SDSS
        field properties for zeropoint, airmass, atmospheric
        extinction, and gain. The noise calculation assumes
        that the number of counts in the converted image is 
        the mean number of counts in the Poisson distribution.
        Thereby, the standard error in that number of counts 
        is the square root of the number of counts in each 
        pixel.
        
        For details on the methods applied here, see:
        http://classic.sdss.org/dr7/algorithms/fluxcal.html
        
        Average quantites obtained from SkyServer SQL form.
        http://skyserver.sdss.org/dr7/en/tools/search/sql.asp
        DR7 Query Form:
        SELECT AVG(airmass_x),AVG(aa_x),AVG(kk_x),AVG(gain_x)
        FROM Field
        '''
        # average sdss photometric field properties (gain is inverse gain)
        airmass  = {'u':1.178, 'g':1.178, 'r':1.177, 'i':1.177, 'z':1.178}
        aa       = {'u':-23.80,'g':-24.44,'r':-24.03,'i':-23.67,'z':-21.98}
        kk       = {'u':0.5082,'g':0.1898,'r':0.1032,'i':0.0612,'z':0.0587}
        gain     = {'u':1.680, 'g':3.850, 'r':4.735, 'i':5.111, 'z':4.622}
        exptime  = 53.907456 # seconds
        # conversion factor from nanomaggies to counts
        counts_per_nanomaggy = exptime*10**(-0.4*(22.5+aa[band]+kk[band]*airmass[band]))
        # image in counts for given field properties
        img_counts = np.clip(img_nanomaggies * counts_per_nanomaggy,a_min=0,a_max=None)
        # poisson noise [adu] computed accounting for gain [e/adu]
        img_counts = np.random.poisson(lam=img_counts*gain[band])/gain[band]
        # convert back to nanomaggies
        img_nanomaggies = img_counts / counts_per_nanomaggy
        
    # add gaussian sky to image
    if common_args['add_false_sky']:
        '''
        Add sky with noise level set by "false_sky_sig" 
        keyword. "false_sky_sig" should be in relative  
        AB magnitudes/arcsec2 units. In other words,
        10**(-0.4*false_sky_sig) gives the sample 
        standard deviation in the sky in linear flux units
        [maggies/arcsec2] around a sky level of zero.
        '''
        # sky sig in AB mag/arcsec2
        false_sky_sig = common_args['false_sky_sig']
        # conversion from mag/arcsec2 to nanomaggies/arcsec2
        false_sky_sig = 10**(0.4*(22.5-false_sky_sig))
        # account for pixel scale in final image
        false_sky_sig *= arcsec_per_pixel**2
        # create false sky image
        sky = false_sky_sig*np.random.randn(*img_nanomaggies.shape)
        # add false sky to image in nanomaggies
        img_nanomaggies += sky
    
    # add image to real sdss sky
    if common_args['add_sdss_sky']:
        '''
        Extract field from galaxy survey database using
        effectively weighted by the number of galaxies in
        each field. For this to work, the desired field
        mask should already have been generated and the
        insertion location selected.
        '''
        import sqlcl
        from astropy.wcs import WCS
        run    = sdss_args['sdss_run']
        rerun  = sdss_args['sdss_rerun']
        camcol = sdss_args['sdss_camcol']
        field  = sdss_args['sdss_field']
        ra     = sdss_args['sdss_ra']
        dec    = sdss_args['sdss_dec']
        exptime = 53.907456 # seconds
        
        # sdss data archive server
        das_url = 'http://das.sdss.org/'
    
        # get and uzip corrected image
        corr_url = das_url+'imaging/{}/{}/corr/{}/'.format(run,rerun,camcol)
        corr_image_name = 'fpC-{:06}-{}{}-{:04}.fit'.format(run,band,camcol,field)
        if not os.access(corr_image_name,0):
            corr_url+='{}.gz'.format(corr_image_name)
            os.system('wget {}'.format(corr_url))
            os.system('gunzip {}'.format(corr_image_name))
        # get wcs mapping
        w = WCS(corr_image_name)
        # determine column and row position in image
        colc,rowc = w.all_world2pix(ra,dec,1,ra_dec_order=True)
        # convert to integers
        colc,rowc = int(np.around(colc)),int(np.around(rowc))
        
        # get field properties from skyServer
        dbcmd = ['SELECT aa_{b},kk_{b},airmass_{b},gain_{b},sky_{b},skysig_{b}'.format(b=band),
                 'FROM Field where run={} AND rerun={}'.format(run,rerun),
                 'AND camcol={} AND field={}'.format(camcol,field)]
        lines = sqlcl.query(' '.join(dbcmd)).readlines()
        # zeropoint, atmospheric extinction, airmass, inverse gain, sky, sky uncertainty
        aa,kk,airmass,gain,sky,skysig = [float(var) for var in lines[1].decode("utf-8").split('\n')[0].split(',')]
        #print(aa,kk,airmass,gain,sky,skysig)
        # convert sky to nanomaggies from maggies/arcsec2
        sky *= (1e9*0.396127**2)
        # convert skysig to nanomaggies from relative sky magnitude errors
        skysig *= sky*np.log(10)/2.5
        # software bias added to corrected images to avoid negative values
        softbias = float(fits.getheader(corr_image_name)['SOFTBIAS'])
        # subtract softbias from corrected image to get image in DN
        corr_image_data = fits.getdata(corr_image_name).astype(float) - softbias # [counts]
        # conversion from nanomaggies to counts
        counts_per_nanomaggy = exptime*10**(-0.4*(22.5+aa+kk*airmass))
        # convert image in counts to nanomaggies with Field properties
        corr_image_data /= counts_per_nanomaggy # [nanomaggies]
        
        if common_args['add_sdss_psf'] and not common_args['add_false_psf']:
            '''
            Grab, reconstruct, and convolve real SDSS PSF image
            with the image in nanomaggies.
            '''
            # get corresponding psf reconstruction image
            psf_url = das_url+'imaging/{}/{}/objcs/{}/'.format(run,rerun,camcol)
            psf_image_name = 'psField-{:06}-{}-{:04}.fit'.format(run,camcol,field)
            if os.access(psf_image_name,0):os.remove(psf_image_name)
            psf_url+=psf_image_name
            os.system('wget {}'.format(psf_url))
            psf_ext = {'u':1,'g':2,'r':3,'i':4,'z':5}
            psfname = 'sdss_psf.fit'
            os.system('{}/Sources/utils/sdss-apps/readAtlasImages-v5_4_11/read_PSF {} {} {} {} {}'.format(realsim_dir,psf_image_name,psf_ext[band],rowc,colc,psfname))
            if os.access(psf_image_name,0): os.remove(psf_image_name)
            # remove softbias from PSF 
            psfdata = fits.getdata(psfname).astype(float)-1000.
            # normalize for convolution with image in nanomaggies
            psfdata /= np.sum(psfdata)
            # convolve with image in nanomaggies
            img_nanomaggies = convolve(img_nanomaggies,psfdata)
            if os.access(psfname,0):os.remove(psfname)
        
        if common_args['add_poisson']:
            '''
            Add Poisson noise to the PSF-convolved image
            with noise level corresponding to the real SDSS
            field properties.
            '''
            # image in counts for given field properties
            img_counts = np.clip(img_nanomaggies * counts_per_nanomaggy,a_min=0,a_max=None)
            # poisson noise [adu] computed accounting for gain [e/adu]
            img_counts = np.random.poisson(lam=img_counts*gain)/gain
            # convert back to nanomaggies
            img_nanomaggies = img_counts / counts_per_nanomaggy
            
        # add real sky pixel by pixel to image in nanomaggies
        corr_ny,corr_nx = corr_image_data.shape
        ny,nx = img_nanomaggies.shape
        for xx in range(nx):
            for yy in range(ny):
                corr_x = int(colc - nx/2 + xx)
                corr_y = int(rowc - ny/2 + yy)
                if corr_x>=0 and corr_x<=corr_nx-1 and corr_y>=0 and corr_y<=corr_ny-1:
                    img_nanomaggies[yy,xx]+=corr_image_data[corr_y,corr_x]
                else:
                    img_nanomaggies[yy,xx]==0.
        if os.access(corr_image_name,0):os.remove(corr_image_name)
        
        # add field info to image header
        warnings.simplefilter('ignore', category=AstropyWarning)            
        header.append(('RUN',run,'SDSS image RUN'),end=True)
        header.append(('RERUN',rerun,'SDSS image RERUN'),end=True)
        header.append(('CAMCOL',camcol,'SDSS image CAMCOL'),end=True)
        header.append(('FIELD',field,'SDSS image FIELD'),end=True)
        header.append(('RA',float(ra),'Cutout centroid RA'),end=True)
        header.append(('DEC',float(dec),'Cutout centroid DEC'),end=True)
        header.append(('COLC',colc,'SDSS image column center'),end=True)
        header.append(('ROWC',rowc,'SDSS image row center'),end=True)
        header.append(('GAIN',gain,'SDSS CCD GAIN'),end=True)
        header.append(('ZERO',aa,'SDSS image zeropoint'),end=True)
        header.append(('EXTC',kk,'SDSS image atm. extinction coefficient'),end=True)
        header.append(('AIRM',airmass,'SDSS image airmass'),end=True)
        header.append(('SKY',sky,'Average sky in full SDSS field [nanomaggies]'),end=True)
        header.append(('SKYSIG',skysig,'Average sky uncertainty per pixel [nanomaggies]'),end=True)
            
    gimage = outputName
    if os.access(gimage,0): os.remove(gimage)
        
#    print('\nAfter Realism:')
#    print('kpc_per_arcsec: {}'.format(kpc_per_arcsec))
#    print('kpc_per_pixel: {}'.format(kpc_per_pixel))
#    print('arcsec_per_pixel: {}'.format(arcsec_per_pixel))
#    m_AB = -2.5*np.log10(np.sum(img_nanomaggies))+22.5
#    print('AB_magnitude: {} at z={}'.format(m_AB,redshift))
#    M_AB = m_AB-5*np.log10(luminosity_distance.value)-25
#    print('AB_Magnitude: {}'.format(M_AB))

    hdu_pri = fits.PrimaryHDU(img_nanomaggies)

    header['REDSHIFT'] = (redshift,'Redshift')
    header.append(('COSMO','FLAT_LCDM','Cosmology'),end=True)
    header.append(('OMEGA_M',cosmo.Om(0),'Matter density'),end=True)
    header.append(('OMEGA_L',cosmo.Ode(0),'Dark energy density'),end=True)
    header.append(('SCALE_1',arcsec_per_pixel,'[arcsec/pixel]'),end=True)
    header.append(('SCALE_2',kpc_per_pixel,'[kpc/pixel]'),end=True)
    header.append(('SCALE_3',kpc_per_arcsec,'[kpc/arcsec]'),end=True)
    header.append(('LUMDIST',cosmo.luminosity_distance(z=redshift).value,'Luminosity Distance [Mpc]'),end=True)
    warnings.simplefilter('ignore', category=AstropyWarning)
    header.extend(zip(common_args.keys(),common_args.values()),unique=True)
    hdu_pri.header = header
    hdu_pri.writeto(gimage)


'''
Script executions start here. This version grabs a corrected image
based on a basis set of galaxies from a database, runs source 
extractor to produce a mask, and selects the location in which to
place the image in the SDSS sky. The final science cutout includes
PSF blurring (real SDSS), SDSS sky from the corrected image, and 
Poisson noise added. The final image is in nanomaggies.
'''

def genSegmap(cutoutName):
    '''Create segmenation image using the sep SExtractor module.'''
    cutoutData = fits.getdata(cutoutName).astype(float)
    # filter kernel
    filter_kernel = np.loadtxt(f'{realsim_dir}/Sources/utils/sdss-cfg/gauss_3.0_7x7.conv',skiprows=2)
    # use std of full image as detection threshold
    guess_rms = np.std(cutoutData)
    # mask all sources above std for background statistics
    mask = ((cutoutData-np.median(cutoutData))>guess_rms)
    # https://github.com/kbarbary/sep/issues/33: convert to float
    # bkg object which includes back() and rms() methods
    bkg = sep.Background(cutoutData, mask=mask, bw=32, bh=32, fw=3, fh=3)
    # run sep.extract() on image
    objCat,segmap = sep.extract(cutoutData-bkg.back(), thresh=1.0, err=bkg.rms(), mask=None, minarea=5,
                             filter_kernel=filter_kernel,filter_type='conv',deblend_nthresh=32,
                             deblend_cont=0.001, clean=True,clean_param=1.0, segmentation_map=True)
    return segmap

def getInjectCoords(segmap):
    '''Use segmentation image to find injection coordinates.
    There is a 10% boundary from the cutout edges which are forbidden.
    A pixels that is both a sky pixel and is inside the boundary is eligible.'''
    # always square cutouts
    nrows = segmap.shape[0]
    ncols = segmap.shape[1]
    # background pixels
    bkgmap = (segmap == 0)
    # pixels within 10% of image size from either side
    bordermap = np.zeros(segmap.shape)
    bordermap[int(0.1*nrows):int(0.9*nrows),int(0.1*ncols):int(0.9*ncols)]=1
    # map of possible injection sites
    segmap = bordermap*bkgmap
    # index of injection site for map
    index = np.random.choice(int(np.sum(segmap)))
    # coordinates of injection site 
    return np.argwhere(segmap)[index]

# get run,rerun,camcol,field,column,row from database data
def rrcf_radec(field_info):
    from astropy.wcs import WCS
    # randomly select from basis set
    index = np.random.randint(low=0,high=len(field_info)-1)
    # define run,rerun,camcol,field of target sky
    run,rerun,camcol,field = field_info[index]
    # sdss data archive server
    das_url = 'http://das.sdss.org/'
    # get and uzip corrected image
    corr_url = das_url+'imaging/{}/{}/corr/{}/'.format(run,rerun,camcol)
    corr_image_name = 'fpC-{:06}-r{}-{:04}.fit'.format(run,camcol,field)
    if not os.access(corr_image_name,0):
        corr_url+='{}.gz'.format(corr_image_name)
        os.system('wget {}'.format(corr_url))
        os.system('gunzip {}'.format(corr_image_name))
    mask_data = genSegmap(corr_image_name)
    rowc,colc = getInjectCoords(mask_data)
    # get wcs mapping
    w = WCS(corr_image_name)
    # determine ra,dec to prevent image registration offsets in each band
    ra,dec = w.all_pix2world(colc,rowc,1,ra_dec_order=True)
    return run,rerun,camcol,field,ra,dec

def make_sdss_args(field_info):
    run,rerun,camcol,field,ra,dec = rrcf_radec(field_info)
    sdss_args = {
                'sdss_run'      : run,   # sdss run
                'sdss_rerun'    : rerun,    # sdss rerun
                'sdss_camcol'   : camcol,     # sdss camcol
                'sdss_field'    : field,   # sdss field
                'sdss_ra'       : ra,  # ra for image centroid
                'sdss_dec'      : dec,   # dec for image centroid
                }
    return sdss_args

# Observational Realism

This repository is for the public release of the statistical observational realism suite described in Bottrell et al (2017ab) and presented publicly in Bottrell et al (2019b). The methods are described in detail in Bottrell et al (2017a) and Bottrell et al (2019b). 

The suite accepts idealized synthetic images in calibrated AB surface brightnesses (Oke & Gunn 1983). Specifically, the input units must be AB mag/arcsec2 (for description of the AB system (see here http://www.sdss3.org/dr8/algorithms/magnitudes.php#nmgy). The standard source for each SDSS band is close to but not exactly the AB source (3631 Jy).

For synthetic images that are not in the rest-frame, the surface brightnesses in each input image bandpass should already be dimmed by (1+z)^-5 for the target redshift. I provide a standalone code which produces idealized synthetic images in AB surface brightnesses for the SDSS bands from SKIRT spectral datacubes (Baes et al 2011, Camps & Baes 2015) covering the optical spectrum.

Images are rebinned to the desired redshift and CCD angular scale. The image must include the physical scale (in pc/pixel) in the FITS header and be identified by the 'CDELT1' keyword. There are then two main options:

(1) The user specifies a point-spread function and sky noise which are incorporated into the image. Examples are provided in which these values are drawn randomly from the distributions of measured sky and PSF in the Sloan Digital Sky Survey. Poisson noise can be added by adopting generic values of photometric calibrations in survey fields (zeropoints, airmass, extinction, CCD gain). 

(2) The images can be inserted into real image fields to incorporate real skies, PSF degradation, and contamination by neighbouring sources in the field of view. The procedure is described in detail in Bottrell et al (2017a) and Bottrell et al (2019b).

The suite may be modified/adapted freely. If you use my suite for your research, I would appreciate a citation to (Bottrell et al (2017a) or the public release in Bottrell et al (2019b).

## Setting up

SourceExtractor must be installed (https://www.astromatic.net/software/sextractor).

An application which reads SDSS PSFs must also be installed. It is included in this package. Go into the Sources/utils/sdss-apps/ directory and do:


rm -rf readAtlasImages-v5_4_11

tar -xzvf readAtlasImages-v5_4_11.tar.gz

cd readAtlasImages-v5_4_11

make clean

make


More info on this package here: https://www.sdss.org/dr12/algorithms/read_psf/

Place the Sources/sqlcl.py file in your Python environment's site-packages directory. Remember to make it executable:

chmod u+x sqlcl.py

Find this script here: http://skyserver.sdss.org/dr7/en/help/download/sqlcl/ for Python 2. Modified to Python 3 for this pipeline by CB. 

We use the Simard et al (2011) quantitative morphologies catalog as a basis for the injection statistics. An SQL version of this catalog can be found at my website here: http://orca.phys.uvic.ca/~cbottrell/share/Realism/sdss_dr7_morph_mybkg_mydeblend_gr.sql. However, we provide a smaller file with only the necessary information with the suite package.

One of the examples in Examples.ipynb (on preparing data from SKIRT datacubes) requires a SKIRT datacube as input to generate idealized gri photometry. You can find one at my website here: http://orca.phys.uvic.ca/~cbottrell/share/Realism/spec_G2G3_e-orbit_1_320_i0_total.fits.



Runs with Python 3 only. Tested with Python 3.6.3, 3.7



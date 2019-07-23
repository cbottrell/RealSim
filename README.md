# `RealSim`

## Statistical observational realism for synthetic images from galaxy simulations

This repository is for the public release of the statistical observational realism suite described in Bottrell et al (2017ab) and presented publicly in Bottrell et al (2019b). The methods are described in detail in Bottrell et al (2017a) and Bottrell et al (2019b). In short, the suite allows one to generate survey-realistic synthetic images of galaxies from hydrodynamical simulations of galaxy formation and evolution. Specifically, the main functionality of this version of `RealSim` inserts "idealized" simulated galaxies into Sloan Digital Sky Survey (SDSS) images in such a way that the statistics of sky brightness, resolution and crowding are matched between simulated galaxies and observed galaxies in the SDSS. 

The suite accepts idealized synthetic images in calibrated AB surface brightnesses (Oke & Gunn 1983). Specifically, the input units must be `AB mag/arcsec2` (for description of the AB system (see here http://www.sdss3.org/dr8/algorithms/magnitudes.php#nmgy). The standard source for each SDSS band is close to but not exactly the AB source (3631 Jy).

For synthetic images that are not in the rest-frame, the surface brightnesses in each input image bandpass should already be dimmed by (1+z)^-5 for the target redshift. I provide a standalone code which produces idealized synthetic images in AB surface brightnesses for the SDSS bands from SKIRT spectral datacubes (Baes et al 2011, Camps & Baes 2015) covering the optical spectrum (see example notebook).

Images are rebinned to the desired redshift and CCD angular scale. The image must include the physical scale (in pc/pixel) in the FITS header and be identified by the `CDELT1` keyword. There are then two main options:

(1) The user specifies a point-spread function and sky noise which are incorporated into the image. Examples are provided in which these values are drawn randomly from the distributions of measured sky and PSF in the Sloan Digital Sky Survey. Poisson noise can be added by adopting generic values of photometric calibrations in survey fields (zeropoints, airmass, extinction, CCD gain). 

(2) The images can be inserted into real image fields to incorporate real skies, PSF degradation, and contamination by neighbouring sources in the field of view. The procedure is described in detail in Bottrell et al (2017a) and Bottrell et al (2019b).

The suite may be modified/adapted freely. If you use my suite for your research, I would appreciate a citation to (Bottrell et al (2017a) and the public release in Bottrell et al (2019b). If you encounter a bug or would like to suggest new functionalities, please contact me and I will promptly get back to you. The `RealSim` methodology can be applied to any existing galaxy imaging survey. I envision multiple versions of this suite -- each tailored to a particular survey's data archive structure, instrumental properties, etc.

## Setting up

### Deblending SDSS corrected images: `Source Extractor`
`Source Extractor` (Bertin & Arnouts 1996) must be installed (https://www.astromatic.net/software/sextractor). I am planning a fix which uses the Pythonized `Source Extractor` module, `sep` (Barbary 2016), for the next big update. This avoids several of external dependencies including configuration files and parameter files as well as the main installation. In preparation for this update, you can do: `pip install sep`. 

### Reading SDSS PSF reconstruction files: `read_psf`
An application which reads SDSS PSFs must also be installed. It is included in this package. Go into the `Sources/utils/sdss-apps/` directory and do:

    rm -rf readAtlasImages-v5_4_11
    tar -xzvf readAtlasImages-v5_4_11.tar.gz
    cd readAtlasImages-v5_4_11
    make clean
    make

More info on this package here: https://www.sdss.org/dr12/algorithms/read_psf/

### Querying the SDSS Data Archive Server (DAS): `sqlcl.py`
Place the `Sources/sqlcl.py` file in your Python environment's site-packages directory. Remember to make it executable:

    chmod u+x sqlcl.py

Find this script here: http://skyserver.sdss.org/dr7/en/help/download/sqlcl/ for `Python 2`. Modified to `Python 3` for this pipeline by CB. 

### Basis catalogue of SDSS galaxies
We use the Simard et al (2011) quantitative morphologies catalogue as a basis for the injection statistics. A `SQL` version of this catalog can be found at my website here: http://orca.phys.uvic.ca/~cbottrell/share/Realism/sdss_dr7_morph_mybkg_mydeblend_gr.sql. However, we provide a smaller file with only the necessary information with the suite package.

### Running the example notebook
One of the examples in `Examples.ipynb` (on preparing data from SKIRT datacubes) requires a SKIRT datacube as input to generate idealized gri photometry. You can find one at my website here: http://orca.phys.uvic.ca/~cbottrell/share/Realism/spec_G2G3_e-orbit_1_320_i0_total.fits. It is 250 MB. It should be added to the `Inputs/Datacubes/` directory to run the `Example.ipynb` notebook.


`RealSim` operates with `Python 3` only. It is tested with versions `3.6` and `3.7`.



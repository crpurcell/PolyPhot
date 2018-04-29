# Tool to perform aperture-photometry on a FITS image

This python script allows the user to display a FITS-format image and
perform aperture-photometry by drawing polygons. Requires some of the
utility scripts in my [pyutils](https://github.com/crpurcell/pyutils)
repository (copy ```util_fits.py```, ```util_phot.py```,
```util_misc.py``` to this directory).

The script was originally created to measure the properties of
extended HII regions detected in
[CORNISH](http://cornish.leeds.ac.uk), but has also been useful for
measuring fluxes in infrared GLIMPSE images. Note that the code
assumes units of Janskys and will divide by the area of a beam, if the
correct headers are present (BMAJ, & BMIN). Details of the formulae
used are given in section 5.3.3 of the [CORNISH catalogue
paper](http://cornish.leeds.ac.uk/public/downloads/Purcell2013_cornish_catalogue.pdf).

The script saves measured values to a text file named for the input
FITS file, except with a '.dat' extension. Multiple measurements can
be saved per image and the script will load previous polygons if
present.

## Usage:
```./polyphot.py G061.4770+00.0891_CORNISH_5GHz.fits```
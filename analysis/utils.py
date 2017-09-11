from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord

import matplotlib.patches as mpatches 

import numpy as np


def SDSS_FITS_filename(RUN, CAMCOL, FIELD, FILTER):

	return "frame-{0}-{1:06d}-{2}-{3:04d}.fits.bz2".format(FILTER, RUN, CAMCOL, FIELD)


def SDSS_spectra_filename(PLATE, MJD, FIBER):
	return 'spec-{0:04d}-{1}-{2:04d}.fits'.format(PLATE, MJD, FIBER)


def open_fits(filename, directory='../data/SDSSfields'):

	hdulist = fits.open('{}/{}'.format(directory, filename))
	return hdulist[0].data, hdulist[0].header


def make_spectra_patch(wcs, ra, dec, radius=8, **kwargs):

	position = SkyCoord(ra, dec, unit='deg').to_pixel(wcs)
	patch = mpatches.Circle(position, radius, fill=False, lw=2, **kwargs)

	return patch


def get_clim(img):

	mean = np.mean(img.ravel())
	std = np.std(img.ravel())
	low, high = (mean-std), (mean+std)*6

	return (low, high)


def open_spectrum(filename, directory):
    # Generate spectrum filename from SDSS_spectra_filename
    # Supply directory which should include galaxy DR7 OBJID
    return fits.open(directory+filename)

def fetch_linelist(spectrum):
    # According to SDSS data structure, the linelist 
    # measurements are in the 4th HDU
    return spectrum[3].data

def fetch_line(linelist, linename):
    # Return all pertinant information about a given line
    return linelist[linelist['LINENAME']==linename]


def SDSS_linedict():
	# Line of interest -- these have to be hardcoded because SDSS
	lines = {'OIII5007': '[O_III] 5007 ',
			 'OIII4959': '[O_III] 4959 ', 
			 'Hb': 'H_beta       ', 
			 'Ha': 'H_alpha      '}

	return lines
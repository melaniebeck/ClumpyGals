from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import astropy.units as u
from astropy.coordinates import SkyCoord


def make_cutout(galdat, image, header, radius='petroR90_r', scale=10):

	gal_position = SkyCoord(galdat['RA'], galdat['DEC'], unit='deg')
	if type(radius) is str: 
		cutout_size = int(scale*galdat[radius]) # scalar --> square cutout
	else:
		cutout_size = radius
	cutout = Cutout2D(image, gal_position, cutout_size, wcs=WCS(header))

	return cutout

def adjust_cutout_hdr(cutout, header):

	newhdr = cutout.wcs.to_header()
	# These params are necessary for pixel scale 
	for key in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
		newhdr[key] = header[key]

	return newhdr


def save_cutout(cutout, header, outfile, **kwargs):

	hdu = fits.PrimaryHDU(cutout.data, header=header)
	hdulist = fits.HDUList(hdu)
	hdulist.writeto("../data/cutouts/"+outfile, **kwargs)


def SDSS_FITS_filename(RUN, CAMCOL, FIELD, FILTER):

	return "frame-{0}-{1:06d}-{2}-{3:04d}.fits.bz2".format(FILTER, RUN, CAMCOL, FIELD)


def open_fits(filename, directory='../data/SDSSfields'):

	hdulist = fits.open('{}/{}'.format(directory, filename))
	return hdulist[0].data, hdulist[0].header


def create_cutout(gal, band='r', **kwargs):
	objid = gal['survey_id']

	rad = 5*gal['petroR90_r']/0.396 # Units: Pixels
	if rad <= 25: rad = 50

	filename = SDSS_FITS_filename(gal['run'], gal['camcol'], 
								  gal['field'], band)
	img, hdr = open_fits(filename)
	cutout = make_cutout(gal, img, hdr, **kwargs)
	chdr = adjust_cutout_hdr(cutout, hdr)

	try: 
		outfile = "{}_{}.fits".format(objid, band)
		save_cutout(cutout, chdr, outfile)
	except:
		pass



from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord

import bz2
import numpy as np
import pdb, glob, os
import argparse


import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches 

import create_cutouts as cut


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


def main(args):
	#maindir = '~/research/clumpygals/data/'

	gals = Table.read("../data/2017-01-30_final_clumpy_sample.csv")
	spectra = Table.read('../data/2017-01-30_object_spectra_lookup_table.csv')
	
	columns = args.columns
	band = args.band
	colors = ['steelblue', 'mediumvioletred', 'coral', 'tomato']

	# Determine the number of rows based on the desired number of columns
	#num = int(round(len(data)/columns)*columns)
	#if num < len(data): num+=columns
	#rows = int(num/columns)*2


	# Loop through chunks of the data set
	for idx in [idx for idx in range(110) if idx%10==0]:
	#for idx in [100]:
		if idx < 100:
			data = gals[idx:idx+10]
		else:
			data = gals[idx:]

		rows = int(round(float(len(data))/2))

		# Start the big ole figure
		fig = plt.figure(figsize=(columns*6,rows*6))
		gs = gridspec.GridSpec(rows, columns, wspace=0., hspace=0.) 


		n = 0
		for i in range(rows):
		#for j in range(columns):
			for j in [0,2]:

				if n < len(data):

					objid = data[n]['survey_id']
					
					infile = "{}_{}.fits".format(objid, band)

					# If the cutout does not already exists, Make It!
					if not os.path.exists("../data/cutouts/"+infile):
						print '../data/cutouts/'+infile+' does not exist!'
						print "Attempting to create cutout from field..."
						cut.create_cutout(data[n], band=band)

					img, hdr = open_fits(infile, directory='../data/cutouts')

					ax = fig.add_subplot(gs[i,j])
					ax.imshow(img, origin='lower', interpolation='nearest',
							  cmap='Greys_r', clim=get_clim(img))
					ax.set_xticks([])
					ax.set_yticks([])

					# Add OBJID, RA, DEC 
					ax.annotate(str(objid), fontsize=16, xy=(0.1, 0.9), 
								color='white', xycoords='axes fraction')
					ax.annotate("{0:.3f}, {1:.3f}".format(data[n]['RA'], 
														  data[n]['DEC']), 
								fontsize=16, xy=(0.1, 0.8), 
								color='white', xycoords='axes fraction')


					"""  ---------   PLOT SPECTRA   ---------  """
					j+=1
					ax2 = plt.subplot(gs[i,j])

					if args.to_plot == 'spectra':
						outstr = 'spectra'

						# Load up all the spectra corresponding to this galaxy
						speclist = glob.glob('../data/SDSSspectra/{0}/*.fits'.format(objid))

						for spec, color in zip(speclist, colors): 
							havespec = True
							try:
								#spectrum = fits.open(filename)
								ss = fits.open(spec)
								spectrum = ss[1]
								lines = ss[3]
							except:
								print "{} not found".format(spec)
								havespec = False


							if havespec:
								# Get the RA, DEC for this spectrum
								ra, dec = ss[0].header['PLUG_RA'], ss[0].header['PLUG_DEC']

								# Create a patch at THAT ^^ position
								wcs = WCS(hdr)
								patch = make_spectra_patch(wcs, ra, dec, 
														   **{'color':color})
								# Add patch to IMG fig
								ax.add_patch(patch)

								flux = spectrum.data['flux']
								lam = 10**spectrum.data['loglam']

								ax2.plot(lam, flux, lw=2., alpha=0.6, color=color)
								ax2.set_xlim(3700, 9500)
								ax2.set_ylim(-50, 1000)

								ax2.annotate("z = {0:.3f}".format(data[n]['z_avg']), 
											fontsize=16, xy=(0.75, 0.9), 
											color='k', xycoords='axes fraction')

					if args.to_plot != 'spectra':
						outstr = 'segmap_'+args.section

						segname = "{}_{}_{}_seg.fits".format(objid, band, 
															 args.section)
						seg = fits.open('../data/cutouts/'+segname)[0].data
						ax2.imshow(seg, origin='lower', clim=get_clim(img))


					ax2.set_xticks([])
					ax2.set_yticks([])
					ax2.set(adjustable='box-forced')
					#"""

				n+=1

		gs.tight_layout(fig)
		plt.savefig('cutouts_{}_{}{}.png'.format(outstr, band, idx))
		#plt.show()


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Display image cutouts along '\
												 'with spectra associated with each galaxy.')
	parser.add_argument('-b', default='r', dest='band',
						help='wavelength band (ugriz) of image you want to view')
	parser.add_argument('-c', default=4, type=int, dest='columns',
						help='Number of columns for resulting image')
	parser.add_argument('-p', default='spectra', dest='to_plot',
						help='Option to plot either the spectra or the seg map')
	parser.add_argument('-s', default='bright', dest='section')

	args = parser.parse_args()
	main(args)
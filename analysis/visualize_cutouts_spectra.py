from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import astropy.units as u
from astropy.coordinates import SkyCoord

import bz2
import numpy as np
import pdb, glob

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches 

import cutout_functions


def SDSS_FITS_filename(RUN, CAMCOL, FIELD, FILTER):

	return "frame-{0}-{1:06d}-{2}-{3:04d}.fits.bz2".format(FILTER, RUN, CAMCOL, FIELD)


def SDSS_spectra_filename(PLATE, MJD, FIBER):
	return 'spec-{0:04d}-{1}-{2:04d}.fits'.format(PLATE, MJD, FIBER)


def open_fits(filename, directory='../data/SDSSfields'):

	hdulist = fits.open('{}/{}'.format(directory, filename))
	return hdulist[0].data, hdulist[0].header


def make_spectra_patch(cutout, ra, dec, radius=8, **kwargs):

	position = SkyCoord(ra, dec, unit='deg').to_pixel(cutout.wcs)
	patch = mpatches.Circle(position, radius, fill=False, lw=2, **kwargs)

	return patch


def get_clim(img):

	mean = np.mean(img.ravel())
	std = np.std(img.ravel())
	low, high = mean-std/4., (mean+std)*2.5

	return (low, high)


def main():
	#maindir = '~/research/clumpygals/data/'

	gals = Table.read("../data/2017-01-30_final_clumpy_sample.csv")
	spectra = Table.read('../data/2017-01-30_object_spectra_lookup_table.csv')
	
	columns = 4
	band = 'r'
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

					try: 
						infile = "{}_{}.fits".format(objid, band)
						img, hdr = open_fits(infile, directory='../data/cutouts')
					except:
						# If the file doesn't exist, create and save the cutout
						pass

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

					#"""
					j+=1
					# Start the spectra figure
					ax2 = plt.subplot(gs[i,j])

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
							patch = make_spectra_patch(cutout, ra, dec, 
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

						ax2.set_xticks([])
						ax2.set_yticks([])
						ax2.set(adjustable='box-forced')
					#"""

				n+=1

		gs.tight_layout(fig)
		plt.savefig('cutouts_spectra{}.png'.format(idx))
		plt.show()


if __name__ == '__main__':
	main()
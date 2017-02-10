import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import glob, pdb
import astropy.io.fits as fits
from astropy.table import Table
import numpy as np
from PIL import Image
import urllib, cStringIO


import SDSSutils.jpegCreator as jc



#function for getting galaxy images
def get_image_from_url(url):

    file = cStringIO.StringIO(urllib.urlopen(url).read())
    img = Image.open(file)
    return img


def plot_img_and_spectra(data, ra='RA', dec='DEC', objid='OBJID', outfile=None):

	rows = int(float(len(data))/2)

	gs = gridspec.GridSpec(rows, 4, wspace=0., hspace=0.5)
	fig = plt.figure(figsize=(len(data), rows*4))

	xy_coords = [(0.02, 0.02),(0.02, 0.97)]


	n = 0
	for i in range(rows):
		for j in [0,2]:
		
			if n < len(data): 

				""" display the galaxy jpeg """

				ax1 = plt.subplot(gs[i,j])

				# get url for the galaxy 
				url = jc.jpeg(ra=data[n][ra], dec=data[n][dec], 
							  height=250, width=250, options='S').url

				# get the image at that url
				img = get_image_from_url(url)

				plt.imshow(img, interpolation='nearest')
				ax1.set_xticks([])
				ax1.set_yticks([])
				ax1.set(adjustable='box-forced')

				ax1.annotate("z = {0:.2}".format(data[n]['z_avg']), fontsize=10, 
							xy=xy_coords[0], color='white', xycoords='axes fraction')

				""" display the spectra associated with this galaxy """

				spectralist = glob.glob('../data/SDSSspectra/{0}/*.fits'.format(data[n][objid]))

				j += 1
				ax2 = plt.subplot(gs[i,j])

				for spectrum in spectralist: 
					spec = fits.open(spectrum)

					flux = spec[1].data['flux']
					lam = 10**spec[1].data['loglam']

					ax2.plot(lam, flux, lw=1., alpha=0.6)
					ax2.set_xlim(3700, 9500)
					#ax2.set_ylim(-50, 1200)

					ax2.set_xticks([])
					ax2.set_yticks([])
					ax2.set(adjustable='box-forced')

			n += 1

	gs.tight_layout(fig)
	plt.savefig('img_and_spec_{}.png'.format(outfile))
	plt.show()



def main():

	dat = Table.read('../data/2017-01-30_final_clumpy_sample.csv')

	for i in [i for i in range(100) if i %10 == 0]:

		if i < 90:
			plot_img_and_spectra(dat[i:i+10], outfile=i, objid='survey_id')
		else:
			plot_img_and_spectra(dat[i:], outfile=i, objid='survey_id')




if __name__ == '__main__':
	main()
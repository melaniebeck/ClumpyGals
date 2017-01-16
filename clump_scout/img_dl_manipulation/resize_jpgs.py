from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import os, sys, pdb
from PIL import Image



def calculate_cropping_box(crop_size, img_size):
	"""
	crop_size will be a percentage of the original image size

	from that percentage, calculate the x,y coordinates of the resulting box
	"""

	x, y = img_size

	x_0, y_0 = img_size[0]/2, img_size[1]/2.

	x1, y1 = int(x/2*(1-crop_size)), int(y/2*(1-crop_size))
	x2, y2 = int(x-x/2*(1-crop_size)), int(y-y/2*(1-crop_size))

	return (x1, y1, x2, y2)




def main():

	clumpy = Table.read('GZ2_no-s82_low-z_low-smooth_clumpzoo.fits')

	# Create a log file to catch any problematic jpgs
	log = open('clumpy_resize.log', 'w')

	'''
	# Create an HTML file to make quick comparison of resizing/cropping
	f = open('clumpy_resize_compare.html','w')
	f.write("""<!doctype html>
<html>
<head>
    <title>Clumpy Galaxies: Crop Only or Resize?</title>
    <link rel="shortcut icon" href="../../img/favicon3.ico" type="image/x-icon">
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
</head>

<body>

<table><tr>
<td width=350 style="text-align:center;">Original</td>
<td width=350 style="text-align:center;">Crop 50%</td>
<td width=350 style="text-align:center;">Crop 70%</td>
<td width=350 style="text-align:center;">Resize 50%</td>
<td width=350 style="text-align:center;">Resize 70%</td>

</tr></table>""")

	'''
    
	for gal in clumpy:
		basename = '{}.jpg'.format(gal['id'])

		# The original jpg will be found in this directory..
		infile = '/data/extragal/willett/gz2/jpg/'+basename

		try:
			img = Image.open(infile)
		except:
			log.write(basename+'\n')

		# Save the image in MY directory
		#img.save('jpg/'+basename)
		#f.write('<img src="{0}" height={s[0]} width={s[1]}>\n'.format('jpg/'+basename, s=img.size))

		if img:
			orig_size = img.size
			crop_size = 0.5

			#for i, crop_size in enumerate(np.arange(0.5, 1.0, 0.1)):
			crop_box = calculate_cropping_box(crop_size, orig_size)

			# Crop original image
			img2 = img.crop(crop_box)
			"""
			#outfile = 'jpg/resized/{0}_crop{1}.jpg'.format(gal['id'], i)
			# Save just the cropped image
			#try: 
			#	img2.save(outfile, "JPEG")
			#	if i==0 or i==3:
			#		f.write('<img src="{0}" height={s[0]} width={s[1]}>\n'.format(outfile, s=img2.size))
			#except:
			#	log.write(outfile+'\n')
			"""

			# Resize the cropped image
			img3 = img2.resize(orig_size, Image.ANTIALIAS)
			outfile = 'jpg/resized/{0}_resized.jpg'.format(gal['id'])

			# Save the resized image -- 
			# I think these will be too blurry!
			try: 
				img3.save(outfile, "JPEG")
				#if (i==0) or (i==3):
				#f.write('<img src="{0}" height={s[0]} width={s[1]}>\n'.format(outfile, s=img3.size))
			except:
				log.write(outfile+'\n')

	log.close()
	f.write("""</body>
</html>""")
	f.close()

if __name__ == '__main__':
	main()

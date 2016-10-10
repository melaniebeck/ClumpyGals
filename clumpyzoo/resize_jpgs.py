from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import os, sys, pdb
from PIL import Image
import glob



def calculate_cropping_box(crop_size, img_size):
	"""
	crop_size will be a percentage of the original image size

	from that percentage, calculate the x,y coordinates of the resulting box
	"""

	x, y = img.size

	x_0, y_0 = img.size[0]/2, img.size[1]/2.

	x1, y1 = int(x/2*(1-crop_size)), int(y/2*(1-crop_size))
	x2, y2 = int(x-x/2*(1-crop_size)), int(y-y/2*(1-crop_size))

	print x1, y1
	print x2, y2

	return (x1, y1, x2, y2)




for infile in sys.argv[1:]:
	img = Image.open(infile)
	orig_size = img.size

	for i, crop_size in enumerate(np.arange(0.5, 1.0, 0.1)):
		crop_box = calculate_cropping_box(crop_size, orig_size)

		img2 = img.crop(crop_box)
		img3 = img2.resize(orig_size, Image.ANTIALIAS)
		outfile = os.path.splitext(infile)[0] + "_crop_resize{}.jpg".format(i)

		try: 
			img3.save(outfile, "JPEG")
		except:
			print "fuuuuuck"
			pdb.set_trace()

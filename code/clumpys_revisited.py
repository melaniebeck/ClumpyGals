from astropy.table import Table, join
import numpy as np
import matplotlib.pyplot as plt
import pdb

#from SDSSutils.jpegs import jpegCreator as jc 

# Import the clumpy galaxies from the coadd stripe82 
clumpy_kw = Table.read('../data/clumpy_stripe82_coadd.fits')

# Import the SDSS portion of the file Mel gave me of ALL 
# clumpy gals with f_clump > 0.5 for GZH
#clumpy_all = Table.read('../data/gzh_sdss.fits')
strictclumpy = Table.read('../data/gzh_sdss_strictclumpy_all_SDSSancillary.fits')

# add U-R color to the table above
uband_corr = strictclumpy['PETROMAG_U']-strictclumpy['EXTINCTION_U']
rband_corr = strictclumpy['PETROMAG_R']-strictclumpy['EXTINCTION_R']

urcolor_corr = uband_corr - rband_corr
urcolor = strictclumpy['PETROMAG_U'] - strictclumpy['PETROMAG_R']

pdb.set_trace()

urllist = []
for objid in clump_all['survey_id']: 
	url = "http://skyserver.sdss.org/dr10/en/tools/explore/obj.aspx?id={0}".format(objid.strip())
	urllist.append(url)

clump_all['location']=urllist
pdb.set_trace()

# remove whitespace of the location column.
urls = [i.strip() for i in clump_all['location']]
clump_all['urls']=urls

# Import the list of urls that I noted as clumpy in my 
# cursory search of the list above
clump_mb = Table.read('blue_clumps_gzh_stripe82.txt', 
                      format='ascii.commented_header')

ind = []
for url in clump_mb['URLS']:
    ind.append(np.where(clump_all['urls']==url)[0])

clump_mb = clump_all[ind]

pdb.set_trace()


# REDO THE ABOVE
# And make manifest script
# RA DEC image name OBJID

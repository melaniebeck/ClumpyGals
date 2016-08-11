#Plots U-R data from fits file 

#Imports necessary packages to read fits file and plot data
from astropy.io import fits
from astropy.table import Table, join
import matplotlib.pyplot as plt
import numpy as np
import pdb
from utils import find_closest




# make region files

def makeregions():
    data = Table.read('clumpy_coadd_sorted.fits')

    blue = data[np.where(data['sample']=='BLUE')]  

    f = open('blue_gals.reg', 'w+')
    f.write('global color=green dashlist=8 3 width=1 font="courier 9 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 include=1 source=1\n')
    f.write('fk5\n')
    for coords in zip(blue['RA'], blue['DEC']):
        f.write('circle(%10.7f,%10.7f,1.2")\n'%(coords[0], coords[1]))
    f.close()
    #pdb.set_trace()

def listfits():
    data = Table.read('clumpy_coadd_sorted.fits')

    data['ra'][0:5]

    blue = data[np.where(data['sample']=='BLUE')]
    
    #,'RUN','RERUN','CAMCOL','FIELD','OBJ','RA','DEC',    
    str1 = 'images/wget das.sdss.org/imaging/' 
    run = blue['RUN']
    camcol = blue['CAMCOL']
    field = blue['FIELD']
    
    for thing in blue['RA']:
        print thing

    
    f = open('blue_gals_fits.sh', 'w+')   
    f.write('ds9s')
    for r, c, d in zip(run, camcol, field):
        if f < 100:
            f.write('images/fpC-%i-g%i-00%i.fit\n'%(r, c, d))
        else:
            f.write('images/fpC-%i-g%i-0%i.fit\n'%(r, c, d))
    
    #pdb.set_trace()
    
def sortsample():
    full = Table.read('clumpy_stripe82_coadd.fits')
    full['sample']='TOSS'
    
    bluesample = Table.read('blue_urls.txt', format='ascii')
    redsample = Table.read('red_urls.txt', format='ascii')
    tosssample = Table.read('toss_urls.txt', format='ascii')
    
    pdb.set_trace()
    
    bluecoords = zip(bluesample['RA'], bluesample['DEC'])
    redcoords = zip(redsample['RA'], redsample['DEC'])
    #tosscoords = zip(toss['RA'], toss['DEC'])
    
    for coord in bluecoords:
        (ra,dec), index = find_closest(coord, zip(full['RA'], full['DEC']))
        full['sample'][index] = 'BLUE'
    for coord in redcoords:
        (ra,dec), index = find_closest(coord, zip(full['RA'], full['DEC']))
        full['sample'][index] = 'RED'
    #for coord in tosscoords:
    #    (ra,dec), index = find_closest(coord, zip(full['RA'], full['DEC']))
    #    full['sample'][index] = 'TOSS'
    full.write('clumpy_stripe82_sorted.fits',overwrite=True)
    
    pdb.set_trace()

def main():
    #makeregions()
    #listfits()
    

    dat = Table.read('clumpy_coadd_sorted.fits')
    mergs = dat[np.where(dat['sample'] == 'RED')]
    ra=mergs['RA']
    dec=mergs['DEC']
    f = open('steve_coords.txt', 'w+')
    i = 1
    for r, d in zip(ra,dec):
        f.write('%i %.5f %.5f\n'%(i, r, d))
        i+=1
    f.close()
    exit()


if __name__ == '__main__':
    main()


'''
str1='http://skyservice.pha.jhu.edu/DR10/ImgCutout/getjpeg.aspx?ra='
str2='&dec='
str3='&scale=0.14&width=512&height=512&opt=&query='

urls = []
for r,d in zip(coadd['RA'], coadd['DEC']):
    url=str1+str(r)+str2+str(d)+str3
    urls.append(url)
    
coadd['urls'] = urls
#--------------------------------------------------------------------------
coadd = Table.read('clumpy_coadd_wurls.fits')

color = coadd['PETROMAG_U']-coadd['PETROMAG_R']
mask = np.ones(len(color), dtype=bool)
mask[np.where(color >= 2.5)] = False

blue = coadd[mask]
red = coadd[~mask]


blues = Table([blue['RA'], blue['DEC'], blue['urls']], names=('ra', 'dec', 'blue urls'))
blues.write('blue_urls.txt', format='ascii.fixed_width')

reds = Table([red['RA'], red['DEC'], red['urls']], names=('ra', 'dec', 'red urls'))
reds.write('red_urls.txt', format='ascii.fixed_width')

blue.write('bluesample_coadd.fits')
red.write('redsample_coadd.fits')
#'''



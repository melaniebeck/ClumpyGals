#Steven Vande Hei
#Plots U-R data from fits file 

#Imports necessary packages to read fits file and plot data
from astropy.io import fits
from pylab import *
import matplotlib as mpl 
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter

hdulist = fits.open('clumpy_stripe82_fullsample.fits')
hdulist1 = fits.open('clumpy_stripe82_coadd.fits')

tbdata1 = hdulist1[1].data
cols = hdulist1[1].columns

tbdata = hdulist[1].data
cols = hdulist[1].columns

#Retrieves U and R values from fits file #1 (clumpy_stripe82_fullsample)
tbdata.field('PETROMAG_U')
tbdata.field('PETROMAG_R')

#Retrieves U and R values from fits file #1 (clumpy_stripe82_coadd)
tbdata1.field('PETROMAG_U')
tbdata1.field('PETROMAG_R')

#Creates new list for all U-R values greater than 2.25
uMinusR=[([x1 - x2 for (x1, x2) in zip(tbdata.field('PETROMAG_U'), tbdata.field('PETROMAG_R'))])]
uMinusR1=[([x1 - x2 for (x1, x2) in zip(tbdata1.field('PETROMAG_U'), tbdata1.field('PETROMAG_R'))])]

#uMinusR=[([x1 - x2 for (x1, x2) in zip(tbdata.field('PETROMAG_U'), tbdata.field('PETROMAG_R')) if (x1 - x2) >= 2.25])]
#uMinusR1=[([x1 - x2 for (x1, x2) in zip(tbdata1.field('PETROMAG_U'), tbdata1.field('PETROMAG_R')) if (x1 - x2) >= 2.25])]


#Modifies tick size 
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2

#Modifies font size
plt.rcParams.update({'font.size': 16})

#Plots frequency histogram uMinusR and uMinusR1
plt.hist(uMinusR, ls='solid', histtype='stepfilled', bins=50, align='right', facecolor='lightsteelblue', edgecolor='black')
plt.hist(uMinusR1, ls='solid', alpha=.5, histtype='stepfilled', bins=50, align='right', facecolor='red', edgecolor='black')
plt.title('Petrosian Flux U - R')
plt.xlabel('U - R')
plt.ylabel('Frequency')
plt.xlim([0,5.5])
plt.show()

#Probability histogram uMinusR and uMinusR1
plt.hist(uMinusR, normed=True, histtype='stepfilled', ls='solid', bins=50, align='right', facecolor='lightsteelblue', edgecolor='black')
plt.hist(uMinusR1, normed=True, alpha=.5, histtype='stepfilled', ls='solid', bins=50, align='right', facecolor='red', edgecolor='black')
plt.title('Petrosian Flux U - R')
plt.xlabel('U - R')
plt.ylabel('Probability')
plt.xlim([0,5.5])
plt.ylim([0, .7])
plt.show()

#Probability accumulation uMinusR and uMinusR1
plt.hist(uMinusR, histtype='stepfilled', normed=True, cumulative=True, ls='solid', bins=50, align='right', facecolor='lightsteelblue', edgecolor='black')
plt.hist(uMinusR1, histtype='stepfilled', alpha=.5, normed=True, cumulative=True, ls='solid', bins=50, align='right', facecolor='red', edgecolor='black')
plt.title('Petrosian Flux U - R')
plt.xlabel('U - R')
plt.ylabel('Accumulation Probability')
plt.xlim([0,5.5])
plt.show()

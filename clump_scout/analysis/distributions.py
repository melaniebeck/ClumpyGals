from __future__ import division

from astropy.table import Table, vstack, join
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
import matplotlib
import pdb
from figure_styles import *

"""
This script can plot the distributions of various galaxy samples for things like
redshift, magnitude, size, etc. It can also create a random subsample that 
has the same distribution as the parent sample. I.E., create a random subsample
of the GZ2 galaxies that has the same distribution in magnitude and size. 
Plots all this shit. Etc.
"""
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()

def gen_bins(data, binsize):
    return np.arange(np.min(data), np.max(data)+binsize, binsize)

def filter_nans(data, columns):
    for column in columns:
        data = data[~np.isnan(data[column])]
    return data
    
def basic_properties(data, subset, filename=None):
    # ---------------- HISTOGRAMS of BASIC PROPERTIES ----------------------
    fig = plt.figure(figsize=(15,5))

    ax1 = fig.add_subplot(131)
    ax1.hist(data['PETROR50_R'], normed=True, histtype='step', lw=2.5, 
             color='grey', bins=gen_bins(data['PETROR50_R'],.5))

    ax1.hist(subset['PETROR50_R'], normed=True, histtype='step', lw=2.5, 
             color='plum', bins=gen_bins(subset['PETROR50_R'],.5))

    ax1.set_xlabel(r'R$_{50,r}$ [arcsec]')
    ax1.set_xlim(0,12)
    
    ax2 = fig.add_subplot(132)
    ax2.hist(data['PETROMAG_R']-data['EXTINCTION_R'], normed=True, 
             histtype='step', lw=2.5, 
             color='k', bins=gen_bins(data['PETROMAG_R'], .1))

    ax2.hist(subset['PETROMAG_R']-subset['EXTINCTION_R'], normed=True, 
             histtype='step', lw=2.5, 
             color='plum', bins=gen_bins(subset['PETROMAG_R'], .1))

    ax2.set_xlabel('m$_r$ [mag]')
    ax2.set_xlim(13,19)
    
    ax4 = fig.add_subplot(133)
    ax4.hist(data['REDSHIFT'], normed=True, histtype='step', lw=2.5, color='k',
             bins=gen_bins(data['REDSHIFT'], .005))

    ax4.hist(subset['REDSHIFT'], normed=True, histtype='step', lw=2.5, color='plum',
             bins=gen_bins(subset['REDSHIFT'], .005))

    xloc = plt.MaxNLocator(6)
    ax4.xaxis.set_major_locator(xloc)
    ax4.set_xlabel('z')
    
    plt.tight_layout()
    plt.savefig('basic_properties_%s.png'%filename)
    plt.show()
    plt.close()


def mass_mag_redshift(data, filename=None):
    data = data[data['MODE']>0] # remove the -1's

    rz_color = data['PETROMAG_G']-data['PETROMAG_Z']
    gi_color = data['PETROMAG_G']-data['PETROMAG_I']

    colors = data['MODE']/13.

    # ----------------- MASS & MAGNITUDE VS REDSHIFT -------------------------
    fig = plt.figure(figsize=(10,10))
    
    ax1 = fig.add_subplot(211)
    plt1 = ax1.scatter(data['REDSHIFT'], data['PETROMAG_R'], 
                       marker='.', c=colors, edgecolor='')
    #ax1.hline(y=17., color='k', ls='--')
    ax1.hlines(y=17.2, xmin=0.02, xmax=0.1, color='k', linestyle='--', lw=2.5)
    ax1.vlines(x=0.1, ymin=10., ymax=17.2, color='k', linestyle='--', lw=2.5)
    ax1.vlines(x=0.02, ymin=10., ymax=17.2, color='k', linestyle='--', lw=2.5)
    ax1.set_xlim(-0.01, .26)
    ax1.set_ylim(10, 18)
    ax1.set_xlabel('z')
    ax1.set_ylabel('m$_r$')
    
    ticks = np.arange(6., 13.)/13.
    #print ticks
    labels = [format(t*13., '.1f') for t in ticks]
    cbax = fig.colorbar(plt1, ticks=ticks)
    cbax.ax.set_ylabel('log M [M/M$_{sun}$]')
    cbax.ax.set_yticklabels(labels)
    
    ax2 = fig.add_subplot(212)
    okay = np.where( (gi_color > 0.) & (gi_color < 2.) )
    
    colors = gi_color[okay]/2.
    
    plt2 = ax2.scatter(data['MODE'][okay], data['REDSHIFT'][okay], marker='.', 
                       c=colors, edgecolor='')
    ax2.hlines(y=0.1, xmin=10., xmax=12.5, color='k', linestyle='--', lw=2.5)
    ax2.hlines(y=0.02, xmin=10., xmax=12.5, color='k', linestyle='--', lw=2.5)
    ax2.vlines(x=10., ymin=0.02, ymax=0.1, color='k', linestyle='--', lw=2.5)
    ax2.set_xlim(6, 12.5)
    ax2.set_ylim(-0.01, .26)
    ax2.set_xlabel('log M [M/M$_{sun}$]')
    ax2.set_ylabel('z')
    
    ticks = np.arange(0., 2.,.5)/2.
    labels = [format(t*2., '.1f') for t in ticks]
    cbax = fig.colorbar(plt2, ticks=ticks)
    cbax.ax.set_ylabel('g - i')
    cbax.ax.set_yticklabels(labels)
    
    plt.tight_layout()
    plt.savefig('mass_mag_z_%s.png'%filename)
    plt.show()
    plt.close()


def select_expert_sample_proportional_GZ2(data):
    subset = data[np.where(data['Expert_label'] != -1)]

    mag_tot = data['PETROMAG_R'] - data['EXTINCTION_R']
    size_tot = data['PETROR50_R']

    mag_sub = subset['PETROMAG_R'] - subset['EXTINCTION_R']
    size_sub = subset['PETROR50_R']

    # -----------------------------------------------------------------------------
    # Divide up the parameter space into bins -- need the bin edges and the 
    # proportion of subjects per bin

    H, xedges, yedges = np.histogram2d(mag_tot, size_tot, normed=True,
                                bins=[gen_bins(mag_tot,.5), gen_bins(size_tot, .5)])

    num_in_subset = 5000.
    running_total = 0.
    h_total = 0.
    for i, x in enumerate(xedges):
        for j, y in enumerate(yedges):
        
            if i < H.shape[0] and j < H.shape[1] and H[i,j]>0.:
                h_total += H[i,j]*.25

                # Select the sample of subjects residing in this bin
                ind = np.where( (mag_tot >= x) & (mag_tot < xedges[i+1]) & 
                                (size_tot >= y) & (size_tot < yedges[j+1]))
                sample = data[ind]

                # How many should we sample from this bin? 
                # (multiply by the bin area to convert from density to probability)
                size = int(round(0.25*H[i,j]*num_in_subset))

                #"""
                if size > 0:
                    rand_ind = rand.randint(0,len(sample), size)
                    running_total += size

                    try:
                        sub_sample = vstack([sub_sample, sample[rand_ind]])
                    except:
                        sub_sample = sample[rand_ind]


    plt.scatter(data['PETROMAG_R']-data['EXTINCTION_R'], data['PETROR50_R'],
                color='k',marker='.', label='Full GZ2')
    plt.scatter(subset['PETROMAG_R']-subset['EXTINCTION_R'], subset['PETROR50_R'],
             color='plum',marker='.', label='Expert')
    plt.scatter(sub_sample['PETROMAG_R']-sub_sample['EXTINCTION_R'], sub_sample['PETROR50_R'],
             color='yellow',marker='.', label='New')

    plt.xlabel('Magnitude')
    plt.ylabel('Size')

    plt.savefig('subsampling_gz2_distributions.png')
    plt.show()

    data = filter_nans(data, ['REDSHIFT'])
    sub_sample = filter_nans(sub_sample, ['REDSHIFT'])
    basic_properties(data, sub_sample, filename='subsample')




###############################################################################
###############################################################################
###############################################################################

gz2 = Table.read('../../GZExpress/SpaceWarps/analysis/GZ2ASSETS_NAIR_MORPH_MAIN.fits')
gz2['OBJID'] = gz2['name']

ancillary = Table.read('../../gzcodez/catalogs_GZ/zoo2AncillaryInfo.fits')

data = join(gz2, ancillary, keys='OBJID')

#smooth = 't01_smooth_or_features_a01_smooth_weighted_fraction'
#feature = 't01_smooth_or_features_a02_features_or_disk_weighted_fraction'
smooth = 't01_smooth_or_features_a01_smooth_debiased'
feature = 't01_smooth_or_features_a02_features_or_disk_debiased'


# --------------------------------------------------------------------------- #
# separate regular SDSS from stripe82
shallow = data[data['stripe82']==0]
deep = data[data['stripe82']==1]

smooth_cut1 = 0.7
smooth_cut2 = 0.8
z_cut1 = 0.06
z_cut2 = 0.1

non_nan = shallow[~np.isnan(shallow['REDSHIFT'])]

print len(non_nan[np.where( (non_nan['REDSHIFT']<=z_cut1) &
							(non_nan[smooth]<=smooth_cut1) ) ])

pdb.set_trace()

clumpy_cut.write('GZ2_no-s82_low-z_low-smooth_clumpzoo.fits')
pdb.set_trace()

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111)

"""
limits = [0., 1., np.min(shallow['PETROMAG_R']), np.max(shallow['PETROMAG_R'])]

counts, extent, aspect = density(shallow[smooth], shallow['PETROMAG_R'], 
                                 limits, numbins=50)
CS = ax.contour(counts.T, extent=extent, levels=[1,10,100,1000], lw=2,
                colors=('plum','mediumorchid','darkorchid','indigo'))#locator=plt.LogLocator()
fmt = ticker.LogFormatterMathtext()
fmt.create_dummy_axis()
manual_locations=[(.8, 11.9),(.81, 12.9), (0.82, 14.1), (0.85, 16.1)]
plt.clabel(CS, CS.levels, fmt=fmt, fontsize=16, manual=manual_locations)

plt.axvline(0.7, ls='--', color='maroon', alpha=.7)
plt.axvline(0.8, ls='--', color='firebrick', alpha=.7)

ax.set_ylim(10., 18.)
ax.set_xlabel('smooth_fraction')
ax.set_ylabel('mag (r band)')
"""

matplotlib.rcParams.update({'font.size': 18, 
							'font.family': 'STIXGeneral', 
							'mathtext.fontset': 'stix'})

# ----------------------------------------- #

#ax = fig.add_subplot(223)

shallow2 = filter_nans(shallow, ['REDSHIFT'])
deep2 = filter_nans(deep, ['REDSHIFT'])

limits = [0., 1., np.min(shallow2['REDSHIFT']), np.max(shallow2['REDSHIFT'])]

counts, extent, aspect = density(shallow2[smooth], shallow2['REDSHIFT'], 
                                 limits, numbins=50)
CS = ax.contour(counts.T, extent=extent, levels=[1,10,100,1000],
                lw=2, colors=('darkgrey','grey', 'black'))
fmt = ticker.LogFormatterMathtext()
fmt.create_dummy_axis()
manual_locations=[(0.7, .12), (0.85, 0.225)]
plt.clabel(CS, CS.levels, fmt=fmt, fontsize=16, manual=manual_locations)


remainder1 = len(shallow[np.where( (shallow['REDSHIFT']<=z_cut1) &
							   	   (shallow[smooth]<=smooth_cut1) ) ])
remainder2 = len(shallow[np.where( (shallow['REDSHIFT']<=z_cut2) &
							   	   (shallow[smooth]<=smooth_cut2) ) ])

ax.plot((smooth_cut1, smooth_cut1), (0., z_cut1), ls='--', lw=3, color='maroon', 
        label='%i remaining'%remainder1)
ax.plot((smooth_cut2, smooth_cut2), (0., z_cut2), ls='-.', lw=3, color='firebrick', 
        label='%i remaining'%remainder2)

ax.plot((0., smooth_cut1), (z_cut1, z_cut1), ls='--', lw=3, color='maroon')
ax.plot((0., smooth_cut2), (z_cut2, z_cut2), ls='-.', lw=3, color='firebrick')

ax.set_xlabel('smooth_fraction')
ax.set_ylabel('redshift')

ax.legend(loc='upper left',frameon=False)



# ----------------------------------------- #
"""
ax = fig.add_subplot(222)

shallow3 = filter_nans(shallow, [smooth])
bins=np.arange(np.min(shallow3[smooth]),np.max(shallow3[smooth])+.05, .05)
ax.hist(shallow3[smooth], histtype='step', color='darkorchid', normed=True, lw=2,
        label='smooth fraction')
remainder1 = np.sum(shallow[smooth]<.7)
remainder2 = np.sum(shallow[smooth]<0.8)
ax.plot((.7, .7), (0., 3.), ls='--', lw=2, color='maroon', 
        label='%i remaining'%remainder1)
ax.plot((.8, .8), (0., 3.), ls='--', lw=2, color='firebrick', 
        label='%i remaining'%remainder2)
ax.legend(loc='best')
"""


"""
ax = fig.add_subplot(222)

limits = [0., 1., np.min(shallow['PETROMAG_R']), np.max(shallow['PETROMAG_R'])]

counts, extent, aspect = density(shallow[feature], shallow['PETROMAG_R'], 
                                 limits, numbins=50)
CS = ax.contour(counts.T, extent=extent, levels=[1,10,100,1000], lw=2,
                colors=('plum','mediumorchid','darkorchid','indigo'))#locator=plt.LogLocator()
fmt = ticker.LogFormatterMathtext()
fmt.create_dummy_axis()
manual_locations=[(.8, 11.9),(.81, 12.9), (0.82, 14.1), (0.85, 16.1)]
plt.clabel(CS, CS.levels, fmt=fmt, fontsize=16, manual=manual_locations)

ax.set_ylim(10., 18.)
ax.set_xlabel('feature_or_disk_fraction')
ax.set_ylabel('mag (r band)')

# ----------------------------------------- #

ax = fig.add_subplot(224)

limits = [0., 1., np.min(shallow2['REDSHIFT']), np.max(shallow2['REDSHIFT'])]

counts, extent, aspect = density(shallow2[feature], shallow2['REDSHIFT'], 
                                 limits, numbins=50)
CS = ax.contour(counts.T, extent=extent, levels=[1,10,100,1000],
                lw=2, colors=('plum','mediumorchid','darkorchid', 'indigo'))
fmt = ticker.LogFormatterMathtext()
fmt.create_dummy_axis()
manual_locations=[(0.7, .12), (0.85, 0.225)]
plt.clabel(CS, CS.levels, fmt=fmt, fontsize=16, manual=manual_locations)

ax.set_xlabel('feature_or_disk_fraction')
ax.set_ylabel('redshift')
"""
plt.tight_layout()
plt.savefig('distribution_of_GZ2_smooths.png')
plt.show()

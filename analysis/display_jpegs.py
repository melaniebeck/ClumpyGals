from __future__ import division

from PIL import Image
import urllib, cStringIO
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pdb
from numpy import random



#function for getting galaxy images
def get_image_from_url(url):

    file = cStringIO.StringIO(urllib.urlopen(url).read())
    img = Image.open(file)
    return img

def show_jpegs(gals, number=None, image='location', criteria=None, 
               randomize=False, order_by=None, title=None, subtitle=None):
    
    if number: 
        num = number
    else:
        num = int(round(len(gals)/5)*5)
        if num < len(gals): num+=5
    
    rows = int(num/5)
 

    if randomize:
        #shuffle which galaxies get displayed   
        galidx = np.arange(len(gals))  
        random.shuffle(galidx)
        gal_list = gals[galidx[:num]]

    else:
        gal_list = gals
    
    if len(gal_list) > num:
        gal_list = gal_list[:num]


    if order_by:
        try:
            gal_list.sort(order_by)
        except:
            print "{0} is not an option for sorting.".format(order_by)


    labels = np.linspace(1,num,num)
    v_align, h_align = ['bottom', 'top'], ['left','left']
    xy_coords = [(0.02, 0.02),(0.02, 0.97)]


    f=plt.figure(figsize=(30,rows*7))
    gs=gridspec.GridSpec(rows,5, wspace=0., hspace=0.)

    n=0
    for i in range (0,rows):
        for j in range(0,5):
            
            ax=plt.subplot(gs[i,j])
            plt.rc('text', usetex=True)

            if n < num:
                gal = gal_list[n] 
            
                plt.imshow(get_image_from_url(gal[image]))
                plt.tick_params(labelbottom='off',labelleft='off')

                # Label each subject with a number in top right corner
                #ax.annotate('%s'%int(labels[n]),fontsize=38,xy=(0.97,0.97),
                #            xycoords='axes fraction', color='white',
                #            verticalalignment='top',horizontalalignment='right')

                try: 
                    for idx, crit in enumerate(criteria):

                        ax.annotate("$\mathrm{%s = %.4f}$"%(crit['label'], 
                                                      gal[crit['colname']]),
                                    fontsize=38,xy=xy_coords[idx], color='white',
                                    xycoords='axes fraction',
                                    verticalalignment=v_align[idx], 
                                    horizontalalignment=h_align[idx])
            
                    ax.set_title(gal['OBJID'], fontsize=25)
                except:
                    pass

                if subtitle:
                    ax.set_title(gal[subtitle], fontsize=25)

                n+=1
                
            else:
                plt.tick_params(labelbottom='off',labelleft='off')

    if title:
        f.text(.5,.92,'\n$f_{%s} > %s$~and~$f_{%s}>%s$'%(display_names[0], values[0], 
                                                         display_names[1], values[1]), 
                                                        fontsize=40, ha='center')
    plt.savefig('clumpyzoo_examples.png')
    plt.show()



def main():

    from astropy.table import Table, join, vstack
    #dat = Table.read('../data/clumpy_coadd_spectra_dr12_gzh.fits')
    #low_z = dat[dat['AVG_Z']<=0.04]

    # construct GZH cat with URLs
    galdat = Table.read('../gz_hubble_main.fits')
    urldat = Table.read('../hubble_main_urls.fits')

    hubble = join(galdat, urldat, keys='zooniverse_id')

    # Isolate Clumpy Gals
    feat = 't01_smooth_or_features_a02_features_or_disk_weighted_fraction'
    clump = 't12_clumpy_a01_yes_weighted_fraction'
    clumpN = 't12_clumpy_total_count'

    clumpy_criteria = ((hubble[feat] > 0.5) & 
                       (hubble[clump] > 0.5) & 
                       (hubble[clumpN] >= 20))

    clumpy = hubble[clumpy_criteria]

    # Select a sample of clumpy gals in redshift bins
    z_bins = np.array([.25, .5, .75, 1., 1.25])

    n = 5
    for i, z in enumerate(z_bins):
        if i < len(z_bins)-1:
            subset = clumpy[ (clumpy['Z_BEST'] >= z) & 
                             (clumpy['Z_BEST'] < z_bins[i+1]) ]

            idx = np.arange(len(subset))
            random.shuffle(idx)
            gals = subset[idx[:n]]

            if i == 0:
                sample = gals
            else:
                sample = vstack([sample, gals])

    # Plot some interesting info for each gal
    stuff = Table.read("""label    colname    value
                 z   Z_BEST      0.04
                 f_{clumpy} t12_clumpy_a01_yes_weighted_fraction 0.5""", format='ascii')
    
    # Create the Image
    show_jpegs(sample, criteria=None, randomize=False, 
               order_by='Z_BEST', image='location')


if __name__ == '__main__':
    main()

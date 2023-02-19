"""
Makes a plot from a FITS file
"""
# Imports
import numpy as np

# The original Python2 version of this document relied on APLPY to load and plot the file, but this was not Python3 compatible.
# Here we have drawn on a script generously provided by George Heald that relies on Astropy instead.

# Updated July 23, 2021 by Mark Richardson and Jayanne English
# Updated Mar 5, 2022 by Gilles Ferrand and Jayanne English; Jayanne English March 5, 2022. 
# Updated Mar 15, 2022 by Mark Richardson and Jayanne English.
# Changes July 8/22: square image changed by Nathan Deg.

#  Changes Dec. 16/22 -- Switched to rectangular trims of image by Nathan Deg
#  Changes Jan. 3/23 -- calculated radius and renamed parameter in "def get_galaxy_range" done by Jayanne Engish
#  Changes Jan. 8/23 -- added figsize as a parameter rather than hard coded. J.E. and Gilles Ferrand. 

import astropy.units as u
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.io import fits
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse, Rectangle

import matplotlib.pylab as pylab

def plot_galaxy(fits_file,RA,DEC,ImgSize,shift,cmap,min_value=None,max_value=None,
                  ticks=None,nsteps=18,label="",coord_frame='fk5',mark_centre=False,show_beam=True,cb_name='',
                  add_tick_ends=True,tick_prec=-2,bkgrd_black=False,title='',TrimSwitch='no_trim', figsize=(8.0,8.0)):

    params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.direction': 'in',
         'ytick.direction': 'in',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
    pylab.rcParams.update(params)

    hdul = fits.open(fits_file)
    hdr = hdul[0].header
    w = WCS(hdr).celestial
    pix_size = np.abs(hdr['CDELT1'])

    # Convert RA to DEG and apply shift
    imagecenterX=15*(RA[0] + RA[1]/60. + RA[2]/3600.) - shift[0]
    imagecenterY=DEC[0] + DEC[1]/60. + DEC[2]/3600 - shift[1]
    imagecenter=[imagecenterX,imagecenterY] 
    w_cut,h_cut=ImageTrim(hdul,w,TrimSwitch, ImgSize, imagecenter,pix_size,coord_frame)
    
    if ticks!=None and add_tick_ends:
        if min_value==None:
            im_min = np.nanmin(h_cut)
            if ticks[0] < im_min:
                min_value = ticks[0]
            else:
                # I want to add a tick to the end of the colour bar. The optional parameter tick_prec sets the precision of this.
                min_tick = (np.floor(im_min/10**tick_prec) + 1)*10**tick_prec
                ticks = [min_tick] + ticks
        if max_value==None:
            im_max = np.nanmax(h_cut)
            if ticks[-1] > im_max:
                max_value = ticks[-1]
            else:
                # I want to add a tick to the end of the colour bar. The optional parameter tick_prec sets the precision of this.
                max_tick = (np.ceil(im_max/10**tick_prec) - 1)*10**tick_prec
                ticks = ticks + [max_tick]

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1,1,1,projection=w_cut)

    cim = ax.imshow(h_cut, cmap=plt.get_cmap(cmap,nsteps), vmin=min_value, vmax=max_value)
    plt.xlabel('RA (J2000)')
    plt.ylabel('Dec (J2000)')
    plt.title(title)

    cbar = fig.colorbar(cim, label=cb_name,ticks=ticks,fraction=0.0467,pad=0.015)

    pix_scale = proj_plane_pixel_scales(w_cut)
    sx, sy = pix_scale[0], pix_scale[1]
    if show_beam:
        try:
            beamx = hdr['BMAJ']/pix_size
            try:
                beamy = hdr['BMIN']/pix_size
                try:
                    beampa = hdr['BPA']
                except:
                    print("No BPA parameter found. Setting position angle to 0 degrees.")
                    beampa = 0.
            except:
                print("No BMIN parameter found. Setting beam to circle.")
                beamy=beamx
                beampa = 0.
            beam = Ellipse((15.,15.), beamx, beamy, angle=beampa,facecolor='black', edgecolor='none', zorder=200)
            ax.add_patch(beam)
        except:
            print("Warning: No beam information found. Beam will not be shown. We suggest setting show_beam=False in plot_galaxy.")

    # Add cross at galaxy centre
    if mark_centre:
        scx,scy = roi.to_pixel(w_cut)
        ax.plot(scx,scy, c='red', marker='+', markersize=12, zorder=300)

    # Make patch in plot black 
    if bkgrd_black:
        ax.set_facecolor((0.0, 0.0, 0.0))

    return fig, ax

def get_galaxy_range(fits_file,RA,DEC,ImgWidth,ImgHeight,shift,coord_frame='fk5'):
    hdul = fits.open(fits_file)
    hdr = hdul[0].header
    w = WCS(hdr).celestial
    pix_size = np.abs(hdr['CDELT1'])

    # Convert RA to DEG and apply shift
    imagecenterX=15*(RA[0] + RA[1]/60. + RA[2]/3600.) - shift[0]
    imagecenterY=DEC[0] + DEC[1]/60. + DEC[2]/3600 - shift[1]

    roi = SkyCoord(imagecenterX, imagecenterY, unit=u.deg, frame=coord_frame)
    pix = skycoord_to_pixel(roi, w)
    RADIUS = np.sqrt(np.square(ImgWidth/2.) + np.square(ImgHeight/2.))
    rsize=np.int(RADIUS/pix_size) + 1

    naxis = len(hdul[0].data.shape)
    if naxis==2: h_cut = hdul[0].data[    int(pix[1]-rsize):int(pix[1]+rsize),int(pix[0]-rsize):int(pix[0]+rsize)] 
    if naxis==3: h_cut = hdul[0].data[0  ,int(pix[1]-rsize):int(pix[1]+rsize),int(pix[0]-size):int(pix[0]+rsize)] 
    if naxis==4: h_cut = hdul[0].data[0,0,int(pix[1]-rsize):int(pix[1]+rsize),int(pix[0]-rsize):int(pix[0]+rsize)] 
    print("Plot range of ", np.nanmin(h_cut),np.nanmax(h_cut))
    return np.nanmin(h_cut),np.nanmax(h_cut)


def ImageTrim(hdul,w,TrimSwitch, TrimLength, imagecenter,pix_size,coord_frame):

    #   Check to make sure we have the correct switch option
    SwitchOptions=['no_trim', 'rectangle']
    if TrimSwitch not in SwitchOptions:
        print("The only valid image triming options are")
        print(TrimSwitch)
        print("Reverting to no trim at all")
        TrimSwitch='no_trim'

    #   First set the center in X and Y coordinates
    imagecenterX=imagecenter[0]
    imagecenterY=imagecenter[1]

    #   Get the center point in pixels rather than RA and DEC
    roi = SkyCoord(imagecenterX, imagecenterY, unit=u.deg, frame=coord_frame)
    pix = skycoord_to_pixel(roi, w)
    
    #   Get the number of axes in the data
    naxis = len(hdul[0].data.shape)
        #   Fix the dimensions of the data to only a 2D slice
    if naxis==2: h_cut = hdul[0].data
    elif naxis==3: h_cut = hdul[0].data[0  ,:,:]
    elif naxis==4: h_cut = hdul[0].data[0,0,:,:]
    
    #   Now get the number of pixels
    npix = [hdul[0].data.shape[-1], hdul[0].data.shape[-2]]
    
    #   No trim option!
    if TrimSwitch=="no_trim":
        h_cut=h_cut
        w_cut = w[0:npix[0]-1,0:npix[1]-1]

    elif TrimSwitch=="rectangle":
        #   Make a 2D set of sizes
        size=np.zeros(2)
        #   Figure out the size of the image in X and Y in pixel units instead of sky-plane units
        size[0]=np.int(TrimLength[0]/pix_size) + 1
        size[1]=np.int(TrimLength[1]/pix_size) + 1
        size=size/2
        #   Make sure that the rectangle does not extend beyond the data image size
        if (int(pix[0]-size[0])<0 or int(pix[0]+size[0])>=npix[0]) or (int(pix[1]-size[1])<0 or int(pix[1]+size[1])>=npix[1]):
            print("WARNING: Requested rectangle is beyond the image size.  Now trimming")
            if int(pix[0]-size[0])<0: size[0] = pix[0]
            if int(pix[0]+size[0])>=npix[0]: size[0] = npix[0]-pix[0]-1
            if int(pix[1]-size[1])<0: size[1] = pix[1]
            if int(pix[1]+size[1])>=npix[1]: size[1] = npix[1]-pix[1]-1
        #   Now we can trim the image to the target size by using pixel indices
        low=pix-size
        high=pix+size
        h_cut = h_cut[int(pix[1]-size[1]):int(pix[1]+size[1]),int(pix[0]-size[0]):int(pix[0]+size[0])]
        #   And do the same to the WCS header portion
        w_cut = w[int(pix[1]-size[1]):int(pix[1]+size[1]),int(pix[0]-size[0]):int(pix[0]+size[0])]
                
        
    return w_cut, h_cut

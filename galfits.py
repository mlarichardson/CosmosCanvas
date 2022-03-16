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

def plot_galaxy(fits_file,RA,DEC,RADIUS,shift,cmap,min_value=None,max_value=None,
                  ticks=None,nsteps=18,label="",coord_frame='fk5',mark_centre=False,show_beam=True,cb_name='',
                  add_tick_ends=True,tick_prec=-2,bkgrd_black=False):
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

    roi = SkyCoord(imagecenterX, imagecenterY, unit=u.deg, frame=coord_frame)
    pix = skycoord_to_pixel(roi, w)
    size=np.int(RADIUS/pix_size) + 1
    npix = [hdul[0].data.shape[-1], hdul[0].data.shape[-2]]
    if (int(pix[0]-size)<0 or int(pix[0]+size)>=npix[0]) or (int(pix[1]-size)<0 or int(pix[1]+size)>=npix[1]): 
        if int(pix[0]-size)<0: size = pix[0]
        if int(pix[0]+size)>=npix[0]: size = npix[0]-pix[0]-1
        if int(pix[1]-size)<0: size = pix[1]
        if int(pix[1]+size)>=npix[1]: size = npix[1]-pix[1]-1
        print("WARNING: Requested radius is larger than data size, trimmed")
    naxis = len(hdul[0].data.shape)
    if naxis==2: h_cut = hdul[0].data[    int(pix[1]-size):int(pix[1]+size),int(pix[0]-size):int(pix[0]+size)] 
    if naxis==3: h_cut = hdul[0].data[0  ,int(pix[1]-size):int(pix[1]+size),int(pix[0]-size):int(pix[0]+size)] 
    if naxis==4: h_cut = hdul[0].data[0,0,int(pix[1]-size):int(pix[1]+size),int(pix[0]-size):int(pix[0]+size)] 
    w_cut = w[int(pix[1]-size):int(pix[1]+size),int(pix[0]-size):int(pix[0]+size)]

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

    fig = plt.figure(figsize=(8.,8.))
    ax = fig.add_subplot(1,1,1,projection=w_cut)

    cim = ax.imshow(h_cut, cmap=plt.get_cmap(cmap,nsteps), vmin=min_value, vmax=max_value)
    plt.xlabel('RA (J2000)')
    plt.ylabel('Dec (J2000)')

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

def get_galaxy_range(fits_file,RA,DEC,RADIUS,shift,coord_frame='fk5'):
    hdul = fits.open(fits_file)
    hdr = hdul[0].header
    w = WCS(hdr).celestial
    pix_size = np.abs(hdr['CDELT1'])

    # Convert RA to DEG and apply shift
    imagecenterX=15*(RA[0] + RA[1]/60. + RA[2]/3600.) - shift[0]
    imagecenterY=DEC[0] + DEC[1]/60. + DEC[2]/3600 - shift[1]

    roi = SkyCoord(imagecenterX, imagecenterY, unit=u.deg, frame=coord_frame)
    pix = skycoord_to_pixel(roi, w)
    size=np.int(RADIUS/pix_size) + 1

    naxis = len(hdul[0].data.shape)
    if naxis==2: h_cut = hdul[0].data[    int(pix[1]-size):int(pix[1]+size),int(pix[0]-size):int(pix[0]+size)] 
    if naxis==3: h_cut = hdul[0].data[0  ,int(pix[1]-size):int(pix[1]+size),int(pix[0]-size):int(pix[0]+size)] 
    if naxis==4: h_cut = hdul[0].data[0,0,int(pix[1]-size):int(pix[1]+size),int(pix[0]-size):int(pix[0]+size)] 
    print("Plot range of ", np.nanmin(h_cut),np.nanmax(h_cut))
    return np.nanmin(h_cut),np.nanmax(h_cut)

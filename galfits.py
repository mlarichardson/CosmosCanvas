"""
Makes a plot from a FITS file
"""
# Imports
import numpy as np

# Loading Gilles backend
# import BackEnd
import aplpy, numpy, math,sys
from astropy.io import fits
from kapteyn import wcs

from matplotlib import pyplot as plt
import matplotlib as mpl

def standard_setup(sp,imagecenterX, imagecenterY, imageradius):
  sp.frame.set_color('black')
  sp.tick_labels.set_font(size='10')
  sp.axis_labels.set_font(size='12')
  sp.set_tick_labels_format(xformat='hh:mm:ss',yformat='dd:mm:ss')
  sp.ticks.set_color('white')
  sp.recenter(imagecenterX, imagecenterY, radius=imageradius)
  #sp.ticks.set_xspacing(0.09)


def colorbar_setup(cb):
  cb.colorbar.set_axis_label_font(size=10)
  cb.colorbar.set_width(0.3)  # arbitrary units, default is 0.2
  cb.colorbar.set_pad(0.07)  # arbitrary units, default is 0.05
  cb.colorbar.set_frame_linewidth(10)
  cb.colorbar.set_font(size='10', \
                      stretch='normal', family='sans-serif', \
                      style='normal', variant='normal')
  cb.colorbar.set_axis_label_pad(2)

def set_color_yr(sp, color):
    '''
    Set the color of the ticks
    '''
    # Major ticks
    for line in sp._ax1.yaxis.get_ticklines():
        line.set_color(color)
    for line in sp._ax2.yaxis.get_ticklines():
        line.set_color(color)

    # Minor ticks
    for line in sp._ax1.yaxis.get_minorticklines():
        line.set_color(color)
    for line in sp._ax2.yaxis.get_minorticklines():
        line.set_color(color)

def plot_galaxy(fits_file,RA,DEC,RADIUS,shift,min_value,max_value,cmap,ticks=None,nsteps=18,label=""):
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'


    # Galaxy and beam observation information:
    #parameters
    alpha_error_ticks = [0.1,0.2,0.3,0.4]
    if ticks == None:
#        alpha_ticks = [-1.4,-1.0,-0.8,-0.6,-0.2,-0.1,0.2]
        alpha_ticks = None
    else:
        alpha_ticks = ticks

    vmin_spix = min_value
    vmax_spix = max_value
    stretch_spix = 'linear'

    # ============================================================= #
    #  INFO ABOVE MUST BE ENTERED FOR EACH GALAXY AND OBSERVATION   #
    # ============================================================= #

    # Nothing below here should be changed
    imagecenterX=15*(RA[0] + RA[1]/60. + RA[2]/3600.) - shift[0]
    imagecenterY=DEC[0] + DEC[1]/60. + DEC[2]/3600 - shift[1]

    imageradius=RADIUS

    fig = plt.figure(figsize=(8.5, 4.3),edgecolor='white',facecolor='white')

    figheight = fig.get_figheight()
    figwidth = fig.get_figwidth()

    subwidth = 1.403
    subheight = 1.2835
    subsize = [subwidth, subheight * (figwidth/figheight) ] # width x height -- this makes it square now y unit is in terms of height r$
    blc = [0.11, 0.04]

    sub_height=2*imageradius # use in the last recenter command
    sub_width=sub_height*(subwidth/subheight) #has to be width divided by 1.429 [0.06 and 0.042]Bigger numbers just zooms out the image$

    # Plots
    fig_plt = aplpy.FITSFigure(fits_file)
    standard_setup(fig_plt,imagecenterX, imagecenterY, imageradius)
    fig_plt.show_colorscale(stretch=stretch_spix, vmin=vmin_spix,vmax=vmax_spix, cmap=plt.get_cmap(cmap,nsteps))

    fig_plt.ticks.set_color('black')
    fig_plt.recenter(x=imagecenterX,y=imagecenterY,width=sub_width,height=sub_height)

    # Colourmaps
    fig_plt.add_colorbar()
    fig_plt.colorbar.set_axis_label_text('alpha')
    colorbar_setup(fig_plt)
    fig_plt.colorbar.set_ticks(alpha_ticks)
    fig_plt.ticks.set_color('black')


    # Labels
    xstart= 0.15
    ystart=0.655

    fig.text(xstart,ystart,label,color='black',size='14',weight='bold', family='serif')

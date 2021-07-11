"""

LIBRARY: specindex
PURPOSE: A Library of tools to generate, compare, test, etc. custom colour maps
    for radio spectral index (alpha) data. Two fixed divergent points are defined
    as alpha_steep <-0.8 and alpha_flat >-0.1. 
AUTHORS: Gilles Ferrand, Mark Richardson, Jayanne English
DATE: Last Edited: 2019-09-13

FUNCTIONS:
    __generate_cmap__
    __draw_path__
    make_cmap_segmented
    __stretch__
    create_cmap_specindex
    create_cmap_velocity
    create_cmap_specindex_error
    __setfig_LCH__
    test_cmap_LCH
    __showme__
    __showme2__
    test_cmap_showme
    test_cmap_showme_compare
"""
# Imports
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from colourspace import convert
from colourspace import gamut
from colourspace import maps

# custom cmap
def __stretch__(p,s1,f1):
    """ If a point has coordinate p  when the steep/flat points are at s0/f0,
        this returns its coordinates when the steep/flat points are at s1/f1.
        This is necessary to preserve the colourmap features.
    """
    s0 = 0.3125
    f0 = 0.75
    dsf = f0 - s0

    if p <= s0:
        return p*s1/s0
    elif p <= f0:
        return ( (p-s0)*f1 + (f0-p)*s1 )/dsf
    else:
        return ( (p-f0) + (1.-p)*f1 )/(1-f0)


def create_cmap_specindex(minp,maxp,steepp=-0.8,flatp=-0.1,name="CC-specindex-default",modes=['clip','crop'],targets=['mpl','png'],png_dir=".",out=False):
    """ Makes a new colour map based on Jayanne English's colourmap
        of yellow - plum, where the orange and dark cyan points
        are fixed to the steep and flat components, while the outer
        regions extend to the min and max values provided.
    """
    # Default:

    # RGB
    # startcolor = '#ffd05f'  # light orange yellow
    # steepcolor = '#d35809'  #  mid orange
    # flatcolor = '#004961'   # dark cyan blue
    # endcolor = '#390358'    # dark magenta-purple i.e. plum
    #
    # minp = -1.3
    # maxp = 0.3
    # steepp = -0.8
    # flatp = -0.1
    #   --> s = (steepp - minp)/(maxp - minp) = 0.3125
    #   --> f = (flatp - minp)/(maxp - minp) = 0.75
    #   --> m = (s+f)/2 = 0.53125
    #
    # The designed colour map defines the step portion in terms of minp - s,
    #  the flat portion in terms of f - maxp, and the intermediate s - p. 
    #  As the user can define all four of these, we need to stretch the 
    #  colourmap as necessary to preserve the features designed into the
    #  colourmap. This is established with the __stretch__ command.

    color_width = maxp - minp

    if steepp < minp or flatp < minp or steepp > maxp or flatp > maxp:
        print("Error: Currently must have minp < steepp < flatp < maxp")
        print("  minp = "), minp
        print("  steepp = "), steepp
        print("  flatp = "), flatp
        print("  maxp = "), maxp
        return None

    s1 = (steepp - minp)/color_width # normalized position of "steep" midpoint
    f1 = (flatp  - minp)/color_width # normalized position of "flat"  midpoint
    m1 = 0.5*(s1+f1)
    
    points = {}
    values = {}

    # Apply stretches to default
    points['L'] = [ 0,   s1, m1, __stretch__(0.6,s1,f1),    f1, __stretch__(0.9,s1,f1),   1]
    values['L'] = [85,   54, 39,     34.3              ,    24,     15.5              ,  15]

    points['C'] = [ 0,   s1, m1, __stretch__(0.6,s1,f1),    f1, __stretch__(0.9,s1,f1),   1]
    values['C'] = [60., 74.4,  0,     7.9               ,  25.1,     46.1              ,  54.4]

    points['H'] = [ 0,   s1, m1, __stretch__(0.6,s1,f1),    f1, __stretch__(0.9,s1,f1),   1]
    values['H'] = [86, 51.7, 72,     200               , 276.2,    302.5              , 320]


    # Apply stretches to default
#    points['L'] = [ 0, s1,  __stretch__(0.5,s1,f1),                      f1, __stretch__(0.9,s1,f1),  1]
#    values['L'] = [85, 54,                      42,                   23.88,                     12, 15]

#    points['C'] = [ 0,   s1, __stretch__(0.527777,s1,f1),                        f1,  1]
#    values['C'] = [60, 74.4,                           0,                     25.14, 60]

#    points['H'] = [ 0,    s1, __stretch__(0.5,s1,f1), __stretch__(0.6,s1,f1),     f1,   1]
#    values['H'] = [86, 51.66,                     72,                    200, 276.24, 320]

    # Apply stretches to default
#    points['L'] = [ 0,   s1, m1, __stretch__(0.6,s1,f1),    f1, __stretch__(0.9,s1,f1),   1]
#    values['L'] = [85,   54, 39,     34.3              ,    23.9,     12.              ,  15]

#    points['C'] = [ 0,   s1, m1, __stretch__(0.6,s1,f1),    f1, __stretch__(0.9,s1,f1),   1]
#    values['C'] = [60., 74.4,  0,     7.9               ,  25.1,     46.1              ,  60.]

#    points['H'] = [ 0,   s1, m1, __stretch__(0.6,s1,f1),    f1, __stretch__(0.9,s1,f1),   1]
#   values['H'] = [86, 51.7, 72,     200               , 276.2,    302.5              , 320]

    if out:
        RGB = maps.make_cmap_segmented(points,values,name=name,modes=modes,targets=targets,png_dir=png_dir,out=out)
        return name, RGB
    else:
        maps.make_cmap_segmented(points,values,name=name,modes=modes,targets=targets,png_dir=png_dir,out=out)
        return name
        

def create_cmap_velocity(minp,maxp,div=0.0, width=0.0):
    """ Makes a color map based on Jayanne English's velocity colourmap
        of yellow - plum, where the orange and dark cyan points
        are fixed to the steep and flat points, while the outer
        regions extend to the min and max values provided.
    """
    # Default:
    color_width = maxp - minp

    d0 = float(div - minp)/float(color_width) # normalized position of midpoint

    points = {}
    values = {}

    # Apply stretches to default
    points['L'] = [ 0,  1]
    values['L'] = [90, 10]

    points['C'] = [ 0, d0, 1]
    values['C'] = [50,  0, 50]

    points['H'] = [ 0, d0-width/2., d0+width/2., 1]

    values['H'] = [30+180, 30+179, 31, 30]

    name_cmap, L,C,H = maps.make_cmap_segmented(points,values,name="blue-red")

    return name_cmap, L, C, H

# testing

def __setfig_LCH__(stack='h'):
    """ Set figures for LCH plots """
    if stack=='h':
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,4))
    if stack=='v':
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(6,9), sharex=True)
        plt.subplots_adjust(hspace=0.1)
    ax1.set_ylabel("L", rotation=0, size=14); ax1.set_xlim(0,1) ; ax1.set_ylim(0,100) ; ax1.set_yticks([0,20,40,60,80,100])
    ax2.set_ylabel("C", rotation=0, size=14); ax2.set_xlim(0,1) ; ax2.set_ylim(0,120) ; ax2.set_yticks([0,30,60,90,120])
    ax3.set_ylabel("H", rotation=0, size=14); ax3.set_xlim(0,1) ; ax3.set_ylim(0,360) ; ax3.set_yticks([0,90,180,270,360])
    return ax1, ax2, ax3

def test_cmap_LCH(L,C,H,Cmax=False,colour=None,width=1.5,stack='h',ax1=None,ax2=None,ax3=None,label=None):
    """ A user test for showing the LCH plots """
    x = np.linspace(0,1,len(L))
    if ax1 == None:
      ax1, ax2, ax3 = __setfig_LCH__(stack=stack)

    ax1.plot(x, L, c=colour, lw=width,label=label)
    if label != None:
      legend = ax1.legend(loc='lower left', shadow=False)

    x = np.linspace(0,1,len(C))
    base_line, = ax2.plot(x, C, c=colour, lw=width)
    if Cmax:
        Cmax = np.zeros(len(C))
        for i in range(len(C)): Cmax[i] = gamut.Cmax_for_LH(L[i],H[i])
        ax2.plot(x, Cmax, c=base_line.get_color(), ls=":", lw=1)
    x = np.linspace(0,1,len(H))
    ax3.plot(x, H%360, c=colour, lw=width)
    return ax1, ax2, ax3

def __showme__(array, cmap="Greys_r", steps=0,minmax_offset=None):
    """ Internal calls to matplotlib to show the colormap image """
    if steps==0: steps = array.shape[1]
    #plt.figure(figsize=(10,4)) # aspect='equal'
    plt.imshow(array, aspect='auto', interpolation='nearest', cmap=plt.cm.get_cmap(cmap,steps))
    plt.gca().set_xticks([],[])
    plt.gca().set_yticks([],[])
    plt.colorbar()
    if minmax_offset == None:
        tiny = (array.max() - array.min())/1000.
    else:
        tiny = minmax_offset
    plt.clim(array.min()-tiny,array.max()+tiny)

def __showme2__(array, cmap_ref="Greys_r", cmap_comp="Greys",steps=0,t1="Ref",t2="Compare",minmax_offset=None):
    """ Internal calls to matplotlib to show a comparison of two colormaps using the colormap image """
    if steps==0: steps = array.shape[1]
    #plt.figure(figsize=(10,4)) # aspect='equal'
    fig = plt.figure(figsize=(18,6))

    ax1 = plt.subplot(121) ; plt.title(t1) ; fig.gca().set_xticks([],[]) ; fig.gca().set_yticks([],[])
    ax2 = plt.subplot(122) ; plt.title(t2) ; fig.gca().set_xticks([],[]) ; fig.gca().set_yticks([],[])


    p1 = ax1.imshow(array, aspect='auto', interpolation='nearest', cmap=plt.cm.get_cmap(cmap_ref,steps))
    if minmax_offset == None:
        tiny = (array.max() - array.min())/1000.
    else:
        tiny = minmax_offset
    p1.set_clim( array.min()-tiny, array.max()+tiny)
    fig.colorbar(p1, ax=ax1)

    p2 = ax2.imshow(array, aspect='auto', interpolation='nearest', cmap=plt.cm.get_cmap(cmap_comp,steps))
    p2.set_clim( array.min()-tiny, array.max()+tiny)
    fig.colorbar(p2, ax=ax2)

def test_cmap_showme(cmap,min_value,max_value,nsteps=18,minmax_offset=None):
    """ A user test to show the colormap image """
    # Make up simple colour comparison

    x = np.arange(-np.pi/2., np.pi/2., 0.01)
    y = np.arange(0,   np.pi, 0.01)
    X,Y = np.meshgrid(x,y)
    Z = min_value + (1 + np.sin(X)*np.sin(Y))/2.*(max_value - min_value)

    __showme__(Z,cmap,nsteps,minmax_offset=minmax_offset)

def test_cmap_showme_compare(cmap1,cmap2,min_value,max_value,nsteps=18,t1="Ref",t2="Compare",minmax_offset=None):
    """ A user test to show a comparison of two colormaps with the colormap images """
    # Make up simple colour comparison

    x = np.arange(-np.pi/2., np.pi/2., 0.01)
    y = np.arange(0,   np.pi, 0.01)
    X,Y = np.meshgrid(x,y)
    Z = min_value + (1 + np.sin(X)*np.sin(Y))/2.*(max_value - min_value)

    __showme2__(Z,cmap1,cmap2,nsteps,t1=t1,t2=t2,minmax_offset=minmax_offset)

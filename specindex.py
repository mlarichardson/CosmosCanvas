"""

LIBRARY: specindex
PURPOSE: A Library of tools to generate, compare, test, etc. custom colour maps
    for radio spectral index (alpha) data. Two fixed divergent points are defined
    as alpha_steep <-0.8 and alpha_flat >-0.1.
AUTHORS: Gilles Ferrand, Mark Richardson, Jayanne English
DATE: Last Edited: 2021-07-11

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


def create_cmap_specindex(min_p,max_p,steep_p=-0.8,flat_p=-0.1,name="CC-specindex-default",modes=['clip','crop'],targets=['mpl','png'],png_dir=".",out=False):
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
    # min_p = -1.3
    # max_p = 0.3
    # steep_p = -0.8
    # flat_p = -0.1
    #   --> s = (steep_p - min_p)/(max_p - min_p) = 0.3125
    #   --> f = (flat_p - min_p)/(max_p - min_p) = 0.75
    #   --> m = (s+f)/2 = 0.53125
    #
    # The designed colour map defines the step portion in terms of min_p - s,
    #  the flat portion in terms of f - max_p, and the intermediate s - p.
    #  As the user can define all four of these, we need to stretch the
    #  colourmap as necessary to preserve the features designed into the
    #  colourmap. This is established with the __stretch__ command.

    color_width = max_p - min_p

    if steep_p < min_p or flat_p < min_p or steep_p > max_p or flat_p > max_p:
        print("Error: Currently must have min_p < steep_p < flat_p < max_p")
        print("  min_p = "), min_p
        print("  steep_p = "), steep_p
        print("  flat_p = "), flat_p
        print("  max_p = "), max_p
        return None

    s1 = (steep_p - min_p)/color_width # normalized position of "steep" midpoint
    f1 = (flat_p  - min_p)/color_width # normalized position of "flat"  midpoint
    m1 = 0.5*(s1+f1)

    # Apply stretches to default
    # Here we set the L,C,H values for the color map. LCH_x correspond to the position along the colourbar (0 at bottom/minimum,
    #  1 at top/maximum). LCH_y is the LCH values at those positions. In between we use a linear interpolation.
    LCH_x_vals = [ 0,    s1, m1, __stretch__(0.6,s1,f1),    f1, __stretch__(0.9,s1,f1),   1]
    LCH_x = {}
    LCH_y = {}

    for coord in ['L','C','H']:
        LCH_x[coord] = np.copy(LCH_x_vals)

    # Each parameter horizontally corresponds to LCH_x_vals positions.
    # Each colour at that position is composed ot the value in L, in C, and in H. So the min_p has the colour 85,60,86.
    # Luminosity ranges from 0 - 100, Chroma ranges from 0 - 100, Hue is degrees on the colour wheel.
    LCH_y['L'] = [85,    54, 39,     34.3              ,    24,     15.5              ,  15]
    LCH_y['C'] = [60., 74.4,  0,     7.9               ,  25.1,     46.1              ,54.4]
    LCH_y['H'] = [86,  51.7, 72,     200               , 276.2,    302.5              , 320]


    if out:
        RGB = maps.make_cmap_segmented(LCH_x,LCH_y,name=name,modes=modes,targets=targets,png_dir=png_dir,out=out)
        return name, RGB
    else:
        maps.make_cmap_segmented(LCH_x,LCH_y,name=name,modes=modes,targets=targets,png_dir=png_dir,out=out)
        return name

def create_cmap_specindex_constantL(L_0=75,C_0=35,H_start=70.,H_dir='left',name="CC-specindex-constL",modes=['clip','crop'],targets=['mpl','png'],png_dir=".",out=False):
    """ Makes a new colour map based on Jayanne English's constant Luminosity/chroma colourmap
        of orange - blue.
    """

    LCH_x = {}
    LCH_y = {}

    for coord in ['L','C','H']:
        LCH_x[coord] = np.arange(0,1.05,0.05)

    # Each parameter horizontally corresponds to LCH_x_vals positions.
    # Each colour at that position is composed ot the value in L, in C, and in H.
    # Luminosity ranges from 0 - 100, Chroma ranges from 0 - 100, Hue is degrees on the colour wheel.
    LCH_y['L'] = [L_0]*len(LCH_x['L'])
    LCH_y['C'] = [C_0]*len(LCH_x['C'])

    if H_dir == 'left':
        H_end = H_start - 180.
    elif H_dir == 'right':
        H_end = H_start + 180.
    else:
        print("Error: H_dir must be 'left' or 'right'")
        return -1

    LCH_y['H'] = [H_start*(1-i) + H_end*i for i in LCH_x['H']]

    if out:
        RGB = maps.make_cmap_segmented(LCH_x,LCH_y,name=name,modes=modes,targets=targets,png_dir=png_dir,out=out)
        return name, RGB
    else:
        maps.make_cmap_segmented(LCH_x,LCH_y,name=name,modes=modes,targets=targets,png_dir=png_dir,out=out)
        return name


def create_cmap_specindex_error(c_mid=0.5,L_ends=72,L_mid=50.,L_min=None,L_max=None,C_max=85.,H_0=70.,H_min=None,H_mid=None,H_max=None,
                                name="CC-specindex-error",modes=['clip','crop'],targets=['mpl','png'],png_dir=".",out=False):
    """ Makes a colour map for uncertainties in spectral index. This is based on Jayanne English's
        error colourmap of light orange and grey, where the pure orange hue indicates the most uncertainty.

        We have added more flexibility for the user, but the default values are setup to produce the desired effect.
    """

    if c_mid < 0. or 1. < c_mid:
        print("Error: Currently must have c_mid outside of [0:1]")
        print("  midp = "), c_mid
        return None

    if L_min==None and L_max==None:
            L_min = L_ends
            L_max = L_ends
    elif L_min==None or L_max==None:
        print("Error: Currently must set both L_min and L_max, or neither and just set L_ends")
        return None

    if H_min==None and H_max==None:
            H_min = H_0
            H_max = H_0
    elif (H_min==None and H_mid==None) or (H_mid==None and H_max==None):
        print("Error: Currently must set either both of H_min and H_max, or H_mid and one of H_mid and H_max, or none of them and just set H_0")
        return None

    C_mid = C_max / 2.
    if H_mid==None:
        H_mid = 0.5*(H_min + H_max)
    elif H_min==None:
        H_min=H_mid
    elif H_max==None:
        H_max=H_mid

    LCH_x_vals = [ 0,  c_mid, 1]
    LCH_x = {}
    LCH_y = {}

    for coord in ['L','C','H']:
        LCH_x[coord] = np.copy(LCH_x_vals)

    # Each parameter horizontally corresponds to LCH_x_vals positions.
    # Each colour at that position is composed ot the value in L, in C, and in H.
    # Luminosity ranges from 0 - 100, Chroma ranges from 0 - 100, Hue is degrees on the colour wheel.
    LCH_y['L'] = [L_min, L_mid, L_max]
    LCH_y['C'] = [0.   , C_mid, C_max]
    LCH_y['H'] = [H_min, H_mid, H_max]


    if out:
        RGB = maps.make_cmap_segmented(LCH_x,LCH_y,name=name,modes=modes,targets=targets,png_dir=png_dir,out=out)
        return name, RGB
    else:
        maps.make_cmap_segmented(LCH_x,LCH_y,name=name,modes=modes,targets=targets,png_dir=png_dir,out=out)
        return name


def create_cmap_velocity(min_p,max_p,div=0.0, width=0.0):
    """ Makes a color map based on Jayanne English's velocity colourmap
        of yellow - plum, where the orange and dark cyan points
        are fixed to the steep and flat points, while the outer
        regions extend to the min and max values provided.
    """
    # Default:
    color_width = max_p - min_p

    d0 = float(div - min_p)/float(color_width) # normalized position of midpoint

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
    cmap = plt.cm.get_cmap(cmap)._resample(steps)
    #plt.figure(figsize=(10,4)) # aspect='equal'
    plt.imshow(array, aspect='auto', interpolation='nearest', cmap=cmap)
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

    cmap_ref = plt.cm.get_cmap(cmap_ref)._resample(steps)
    cmap_comp = plt.cm.get_cmap(cmap_comp)._resample(steps)
    #plt.figure(figsize=(10,4)) # aspect='equal'
    fig = plt.figure(figsize=(18,6))

    ax1 = plt.subplot(121) ; plt.title(t1) ; fig.gca().set_xticks([],[]) ; fig.gca().set_yticks([],[])
    ax2 = plt.subplot(122) ; plt.title(t2) ; fig.gca().set_xticks([],[]) ; fig.gca().set_yticks([],[])


    p1 = ax1.imshow(array, aspect='auto', interpolation='nearest', cmap=cmap_ref)
    if minmax_offset == None:
        tiny = (array.max() - array.min())/1000.
    else:
        tiny = minmax_offset
    p1.set_clim( array.min()-tiny, array.max()+tiny)
    fig.colorbar(p1, ax=ax1)

    p2 = ax2.imshow(array, aspect='auto', interpolation='nearest', cmap=cmap_comp)
    p2.set_clim( array.min()-tiny, array.max()+tiny)
    fig.colorbar(p2, ax=ax2)

def test_cmap_showme(cmap,min_value=0,max_value=1,nsteps=18,minmax_offset=None):
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

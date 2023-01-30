"""

LIBRARY: specindex
PURPOSE: A Library of tools to generate, compare, test, etc. custom colour maps
    for radio spectral index (alpha) data. Two fixed divergent points are defined
    as alpha_steep <-0.8 and alpha_flat >-0.1.
AUTHORS: Gilles Ferrand, Mark Richardson, Jayanne English
DATE: Last Edited: 2023-01-30 GF.

FUNCTIONS:
    __stretch__
    create_cmap_specindex
    create_cmap_specindex_constantL
    create_cmap_specindex_error
"""
# Imports
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from colourspace import convert
from colourspace import gamut
from colourspace import maps
from sys import exit

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


def create_cmap_specindex(min_p,max_p,steep_p=-0.8,flat_p=-0.1,name="CC-specindex-default",mode='clip',targets=['mpl'],mpl_reg=True,out=False):
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

    # Check mode is a string or single element list, otherwise issue a warning
    if isinstance(mode,str):
        modes = [mode]
    elif isinstance(mode,list):
        modes=mode
        if len(mode)>1:
            print("Warning: ColourConvas tutorials only address a single mode of colourmap from colourspace (either 'clip' or 'crop').")
            print("Warning: By providing both, the colour map names will match the 'name' argument with suffix '_clip' and '_crop'.")
            print("Warning: Please ensure that you wish to use colourspace and CosmosCanvas in this way. Expected 'mode' is a string.")
    else:
        print("Error. Expected 'mode' to be a string. 'mode' can also be a list. 'mode' has value and type:", mode, type(mode))
        exit(-1)

    RGB = maps.make_cmap_segmented(LCH_x,LCH_y,name=name,modes=modes,targets=targets,mpl_reg=mpl_reg,out=out)
    if out: return RGB

def create_cmap_specindex_constantL(L_0=75,C_0=35,H_start=70.,H_dir='left',name="CC-specindex-constL",mode='clip',targets=['mpl'],mpl_reg=True,out=False):
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

    # Check mode is a string or single element list, otherwise issue a warning
    if isinstance(mode,str):
        modes = [mode]
    elif isinstance(mode,list):
        modes=mode
        if len(mode)>1:
            print("Warning: ColourCanvas tutorials only address a single mode of colourmap from colourspace (either 'clip' or 'crop').")
            print("Warning: By providing both, the colour map names will match the 'name' argument with suffix '_clip' and '_crop'.")
            print("Warning: Please ensure that you wish to use colourspace and CosmosCanvas in this way. Expected 'mode' is a string.")
    else:
        print("Error. Expected 'mode' to be a string. 'mode' can also be a list. 'mode' has value and type:", mode, type(mode))
        exit(-1)

    RGB = maps.make_cmap_segmented(LCH_x,LCH_y,name=name,modes=modes,targets=targets,mpl_reg=mpl_reg,out=out)
    if out: return RGB


def create_cmap_specindex_error(c_mid=0.5,L_ends=72,L_mid=50.,L_min=None,L_max=None,C_max=85.,H_0=70.,H_min=None,H_mid=None,H_max=None, name="CC-specindex-error",mode='clip',targets=['mpl'],mpl_reg=True,out=False):
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

    # Check mode is a string or single element list, otherwise issue a warning
    if isinstance(mode,str):
        modes = [mode]
    elif isinstance(mode,list):
        modes=mode
        if len(mode)>1:
            print("Warning: ColourConvas tutorials only address a single mode of colourmap from colourspace (either 'clip' or 'crop').")
            print("Warning: By providing both, the colour map names will match the 'name' argument with suffix '_clip' and '_crop'.")
            print("Warning: Please ensure that you wish to use colourspace and CosmosCanvas in this way. Expected 'mode' is a string.")
    else:
        print("Error. Expected 'mode' to be a string. 'mode' can also be a list. 'mode' has value and type:", mode, type(mode))
        exit(-1)
        
    RGB = maps.make_cmap_segmented(LCH_x,LCH_y,name=name,modes=modes,targets=targets,mpl_reg=mpl_reg,out=out)
    if out: return RGB


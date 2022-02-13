"""

LIBRARY: specindex
PURPOSE: A Library of tools to generate, compare, test, etc. custom colour maps
    for radio spectral index (alpha) data. Two fixed divergent LCH_x are defined
    as alpha_steep <-0.8 and alpha_flat >-0.1.
AUTHORS: Gilles Ferrand, Mark Richardson, Jayanne English
DATE: Last Edited: 2022-01-07

FUNCTIONS:
    __stretch__
    create_cmap_specindex
    create_cmap_specindex_constantL
    create_cmap_specindex_error
    create_cmap_velocity
"""
# Imports
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
from colourspace import convert
from colourspace import gamut
from colourspace import maps
from sys import exit
try:
    from pathlib import Path
except:
    print("Warning: pathlib module not available. Will not check that png_dir exists.")

# custom vel cmap
def create_cmap_velocity(min_p,max_p,div=0.0, width=0.001,Lval_max=90.,Lpoint_1=0.33,Lval_1=61,Lval_2=None,Lval_mid=None,
      Lval_min=10.,Cval_max=50,Cval_1=None,Cval_2=0.4,Hval_L=210.,Hval_R=30.,Hval_1=210.,Hval_2=209,Hval_3=31,Hval_4=30., name="blue-red",
      mode='clip',targets=['mpl','png'],mpl_reg=True,png_dir="./cmaps",out=False):
    """ Makes a color map based on Jayanne English's velocity colourmap.
        This allows a lot of freedom in how the L,C,H varies.
        L:
            - L transitions linearly from (0,Lval_max) --> (Lpoint_1,Lval_1) --> (div - width/2,Lval_2) --> (div,Lval_mid) (symmetric linear about div)
        C:
            - C transitions linearly from (0,Cval_max) --> (Lpoint_1,Cval_1) --> (div - width/2,Cval_2) --> (div,0.) (symmetric reflected about div)
        H:
            - H transitions linearly from (0,Hval_L) --> (Lpoint_1,Hval_1) --> (div - width/2,Hval_2) --> (div+width/2,Hval_3) --> (1 - Lpoint_1,Hval_4) --> (1,Hval_R)

    """
    # Default:
    color_width = max_p - min_p

    d0 = float(div - min_p)/float(color_width) # normalized position of midpoint
    p2 = d0 - width/2.
    p3 = d0 + width/2.

    if Lval_mid == None:
        Lval_mid = (Lval_max + Lval_min)/2.
    if Lval_2 == None:
        Lval_2 = ((p2-Lpoint_1)*Lval_mid + (d0-p2)*Lval_1)/(d0-Lpoint_1)
        (d0-width/2.)*(Lval_1-Lval_max)/Lpoint_1 + Lval_max
    if Cval_1 == None:
        Cval_1 = Cval_max
    Lval_3 = Lval_min + (Lval_max - Lval_2)
    Lval_4 = Lval_min + (Lval_max - Lval_1)
    Hval_mid = (Hval_2 + Hval_3)/2.

    LCH_x = {}
    LCH_y = {}

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

    # Apply stretches to default
    LCH_x['L'] = np.array([       0.,  Lpoint_1, d0-width/2.,       d0, d0+width/2., 1-Lpoint_1,        1.])
    LCH_y['L'] = [Lval_max,    Lval_1,      Lval_2, Lval_mid,      Lval_3,     Lval_4, Lval_min]

    LCH_x['C'] = np.copy(LCH_x['L'])
    LCH_y['C'] = [Cval_max,     Cval_1,     Cval_2,       0.,      Cval_2,      Cval_1, Cval_max]

    LCH_x['H'] = np.copy(LCH_x['L'])
    LCH_y['H'] = [Hval_L,       Hval_1,     Hval_2,  Hval_mid,     Hval_3,      Hval_4,    Hval_R]

    RGB = maps.make_cmap_segmented(LCH_x,LCH_y,name=name,modes=modes,targets=targets,mpl_reg=mpl_reg,png_dir=png_dir,out=out)
    if out: return RGB

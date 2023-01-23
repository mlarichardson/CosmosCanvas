"""

LIBRARY: velmap
PURPOSE: A Library of tools to generate, compare, test, etc. custom colour maps
        for velocity map (moment 0) data. Several divergent points can be set, 
        although the default is a map works creates a colour legend that works with human perception.
        It consists of 2 complementarly colours (cyan blue == blueshift; red-orange == redshift) 
        and luminosity that monotonically increases from redshift to blueshift (low values are lightest). 
AUTHORS: Gilles Ferrand, Mark Richardson, Jayanne English
DATE: Last Edited: 2022-12-16 JE

FUNCTIONS:
    __stretch__
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

# custom Single Divergent Point velocity cmap
# JE: added to def Hval_mid=None and Cval_mid=0.0

def create_cmap_velocity(min_p,max_p,div=0.0, width=0.05,Lval_max=90.,Lpoint_1=0.33,Lval_1=61,Lval_2=None,Lval_mid=None,
 Lval_min=10.,Cval_max=50,Cval_mid=None,Cval_1=None,Cval_2=0.4,Hval_mid=None,Hval_L=210.,Hval_R=30.,Hval_1=210.,Hval_2=209,Hval_3=31,Hval_4=30., Lval_3=None, Lval_4=None, name="blue-red", 
      mode='clip',targets=['mpl','png'],mpl_reg=True,png_dir="./cmaps",out=False):
    """ Makes a color map based on Jayanne English's velocity colourmap.
        This allows a lot of freedom in how the L,C,H varies.
        L:
            - L transitions linearly from (0,Lval_max) --> (Lpoint_1,Lval_1) --> (div - width/2,Lval_2) --> (div,Lval_mid) (symmetric linear about div)
        C:
            - C transitions linearly from (0,Cval_max) --> (Lpoint_1,Cval_1) --> (div - width/2,Cval_2) --> (div,0.) (symmetric reflected about div)
        H:
            - H transitions linearly from (0,Hval_L) --> (Lpoint_1,Hval_1) --> (div - width/2,Hval_2) --> (div+width/2,Hval_3) --> (1 - Lpoint_1,Hval_4) --> (1,Hval_R)

        Width: The 'width' parameter is the region where chroma is zero (i.e. grey) centred on div and lies between Lpoint_2 and Lpoint_3. Note Lval_mid is at the middle of this region and therefore if the width is too small the Lval_mid may not appear to reach the set luminosity value. NB: for large widths an optical illusion occurs that causes one to see a complementary colour where there actually exists grey. For example, if the grey is near orange you may see cyan. 

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
    if Lval_3 == None:
        Lval_3 = Lval_min + (Lval_max - Lval_2)
    if Lval_4 == None:
        Lval_4 = Lval_min + (Lval_max - Lval_1)
        print("Cval_mid is ")
        print(Cval_mid)
    if Cval_mid == None:
        Cval_mid = 0.0
    if Hval_mid == None:    
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
    #LCH_y['C'] = [Cval_max,     Cval_1,     Cval_2,       0.,      Cval_2,      Cval_1, Cval_max] 
    #JE: changed for customizing. Default should be ZERO. 
    LCH_y['C'] = [Cval_max,     Cval_1,     Cval_2,       Cval_mid,      Cval_2,      Cval_1, Cval_max]

    LCH_x['H'] = np.copy(LCH_x['L'])
    LCH_y['H'] = [Hval_L,       Hval_1,     Hval_2,  Hval_mid,     Hval_3,      Hval_4,    Hval_R]

    RGB = maps.make_cmap_segmented(LCH_x,LCH_y,name=name,modes=modes,targets=targets,mpl_reg=mpl_reg,png_dir=png_dir,out=out)
    
    print('Value of Cval_mid is ')
    print(LCH_y['C'])
    
    if out: return RGB
    
    
    
def create_cmap_doubleVelocity(minvalue,maxvalue,div=0.0, Cval_max=35, name="CC-vmap-DoubleLum1-default"):
        
    # Parameters for Hues:
    # Hval_L,       Hval_1,     Hval_2,  Hval_mid,     Hval_3,      Hval_4,    Hval_R
    # For lowest data value (blueshift) at "L"=left through to highest redshift at "R"=right.

    # Similarly change luminosities with the parameters
    # LCH_y['L'] = [Lval_max,    Lval_1,      Lval_2, Lval_mid,      Lval_3,     Lval_4, Lval_min]
    
    # Select Hues, using degrees on the colour wheel: 
    Hval_L=190. # Left for lowest data value.  Turquoise.
    Hval_R=10.  # Right for highest data value. Rose.
    Hval_1=210. # cyans
    Hval_2=230.
    Hval_3=40.  # Set this point, which is after grey midpoint, to higher (redshift) data values. 
                # While hue 50 is the complement to 230 but looks too brown.
    Hval_4=30.

    # Adjust luminosity: 
    Lval_1=61. # This is the default value in .py file.
    #Lval_2=55. #setting a value here does not permit a change in Lval_mid=85. 
    Lval_2=None
    Lval_mid=55.
    Lval_3=40 
    Lval_4=30
    
    VMap2= create_cmap_velocity(minvalue,maxvalue,Cval_max=Cval_max,div=div,name=name,
                                Hval_L=Hval_L,Hval_R=Hval_R,Hval_1=Hval_1,Hval_2=Hval_2,Hval_3=Hval_3,
                                Hval_4=Hval_4, Lval_mid=Lval_mid,Lval_1=Lval_1,Lval_2=Lval_2,Lval_3=Lval_3, Lval_4=Lval_4)

    return name
# -*- coding: utf-8 -*-
"""
Illumination Correction

Notes:
    -  corrects image via (img - bkg) / (flt - bkg) 

Reference: 
    - Fundamentals of Light Microscopy and Electronic Imaging, p 421
    
@author: ckirst
"""


import sys

import iDISCO.IO.IO as io

from iDISCO.Utils.Timer import Timer
from iDISCO.Utils.ParameterTools import writeParameter

from iDISCO.Visualization.Plot import plotTiling

def correctIllumination(img, flatfield, background = None, subStack = None, verbose = False, illuminationFile = None, out = sys.stdout, **parameter):
    """Correct Illumination"""
    timer = Timer();
    if flatfield is None or isinstance(flatfield, str):
        fld = flatfield;
    else:
        fld = "image of size %s" % str(flatfield.shape);
    
    if background is None or isinstance(background, str):
        bkg = background;
    else:
        bkg = "image of size %s" % str(background.shape);
    
    writeParameter(out = out, head = 'Illumination', flatfield = fld,  background = bkg);
    
    #read data    
    if flatfield is None:
        return img;
    flatfield = io.readData(flatfield);
    background = io.readData(background);
    
    if flatfield.shape != img[:,:,0].shape:
        raise RuntimeError("correctIllumination: flatfield does not match image size: %s vs %s" % (flatfield.shape,  img[:,:,0].shape));
    
    # illumination correction in each slice
    if background is None:
        for z in range(img.shape[2]):
            img[:,:,z] = img[:,:,z] / flatfield;
    else:
        if background.shape != flatfield.shape:
            raise RuntimeError("correctIllumination: background does not match image size: %s vs %s" % (background.shape,  img[:,:,0].shape));        

        fld = (flatfield - background);
        for z in range(img.shape[2]):
            img[:,:,z] = (img[:,:,z] - background) / fld;
    
    #write result for inspection
    if not illuminationFile is None:
        if not subStack is None:
            ii = subStack["zSubStackCenterIndices"][0];
            ee = subStack["zSubStackCenterIndices"][1];
            si = subStack["zCenterIndices"][0];
        else:
            si = 0;
            ii = 0;
            ee = -1;
        
        io.writeData(illuminationFile, img[:,:,ii:ee], startIndex = si );     
     
    #plot result for inspection
    if verbose:
        plotTiling(img);

    out.write(timer.elapsedTime(head = 'Illumination') + '\n');    
    return img 
    













   
## stuff
#   
#from iDISCO.Analysis.Tools.Extrapolate import extrap1d 
#from scipy.interpolate import interp1d
#def flatfieldFromIntensityProfile(profile, size, axis = 1):
#    """Generates a flatfield image from a 1d intensity profile along major axis of illumination distortion"""
#    
#    size = io.dataSize(size);
#    size = size[0:2];
#    
#    # profile is array of x,intensity
#    pf = extrap1d(profile[:,0], profile[:,1], interpolation = 'cubic', exterpolation = 'constant');
#    pf = pf(numpy.array(range(0, size[axis])));
#    
#    if axis == 1:
#        
#        else:
#    
#    
#    
#>>>
#
#>>> x = np.linspace(0, 10, num=11, endpoint=True)
#>>> y = np.cos(-x**2/9.0)
#>>> f = interp1d(x, y)
#>>> f2 = interp1d(x, y, kind='cubic')
#    
#    
#xnew = np.linspace(0, 10, num=41, endpoint=True)
#>>> import matplotlib.pyplot as plt
#>>> plt.plot(x, y, 'o', xnew, f(xnew), '-', xnew, f2(xnew), '--')
#>>> plt.legend(['data', 'linear', 'cubic'], loc='best')
#>>> plt.show()    
#    
#    
#def backgroundFromOpening()

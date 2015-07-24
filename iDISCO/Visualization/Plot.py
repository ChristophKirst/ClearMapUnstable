# -*- coding: utf-8 -*-
"""
Plotting routines for overlaying lables, tilings, and sectioning of 3d data sets

Created on Fri Jun  5 00:04:29 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];


import math
import numpy
import matplotlib as mpl
import matplotlib.pyplot as plt


def plotTiling(image, tiling = "automatic", maxtiles = 20): 
    """Plot 3d image as tiles"""
    
    dim = len(image.shape);
     
    if dim < 3 or dim > 4:
        raise StandardError('plotTiling: image dimension must be 3  or 4');
    
    if len(image.shape) == 3: 
        if image.shape[2] == 3:  # 2d color image
            ntiles = 1;
            cmap = None;
            image = image.reshape((image.shape[0], image.shape[1], 1, 3));
        else:                   # 3d gray image
            ntiles = image.shape[2];
            cmap = plt.cm.gray;
            image = image.reshape(image.shape + (1,));
    else:
        ntiles = image.shape[2]; # 3d color = 4d
        cmap = None;
    
    if ntiles > maxtiles:
        print "plotTiling: number of tiles %d very big! Clipping at %d!" % (ntiles, maxtiles);
        ntiles = maxtiles;
    
    
    if tiling == "automatic":
        nx = math.floor(math.sqrt(ntiles));
        ny = int(math.ceil(ntiles / nx));
        nx = int(nx);
    else:
        nx = int(tiling[0]);
        ny = int(tiling[1]);     
    
    #print image.shape
        
    fig, axarr = plt.subplots(nx, ny, sharex = True, sharey = True);
    axarr = numpy.array(axarr);
    axarr = axarr.flatten();
    
    imin = image.min();
    imax = image.max();
    
    for i in range(0, ntiles): 
        a = axarr[i];
        imgpl = image[:,:,i,:];
        if imgpl.shape[2] == 1:
            imgpl = imgpl.reshape((imgpl.shape[0], imgpl.shape[1]));       
        a.imshow(imgpl, interpolation='none', cmap = cmap, vmin = imin, vmax = imax);

    fig.canvas.manager.window.activateWindow()
    fig.canvas.manager.window.raise_()
    
    return fig;
 
def overlayLabel(image, label, alpha = False, lcmap = 'jet'):
    """Overlay a gray scale image with colored labeled image"""
    
    lmax = label.max();
    
    if lmax <= 1:
        carray = numpy.array([[1,0,0,1]]);
    else:
        cm = mpl.cm.get_cmap(lcmap);
        cNorm  = mpl.colors.Normalize(vmin=1, vmax = int(lmax));
        carray = mpl.cm.ScalarMappable(norm=cNorm, cmap=cm);
        carray = carray.to_rgba(numpy.arange(1, int(lmax + 1)));

    if alpha == False:
        carray = numpy.concatenate(([[0,0,0,1]], carray), axis = 0);
    else:
        carray = numpy.concatenate(([[1,1,1,1]], carray), axis = 0);
        
    cm = mpl.colors.ListedColormap(carray);
    carray = cm(label);
    carray = carray.take([0,1,2], axis = -1);

    if alpha == False:
        cimage = (label == 0) * image;
        cimage = numpy.repeat(cimage, 3);
        cimage = cimage.reshape(image.shape + (3,));        
        cimage += carray;
    else:
        cimage = numpy.repeat(image, 3);
        cimage = cimage.reshape(image.shape + (3,));    
        cimage *= carray;

    return cimage;
    
    
def plotOverlayLabel(image, label, alpha = False, tiling = "automatic", maxtiles = 20):
    """Plot gray scale image overlayed with labeled image"""
    ov = overlayLabel(image, label, alpha = alpha);
    plotTiling(ov, tiling = tiling);


    
def plotOverlayPoints(image, points, alpha = False, tiling = "automatic", maxtiles = 20):
    
    dsize = list(image.shape);
    if dsize[2] > maxtiles:
        dsize[2] = maxtiles;
    
    imgc = numpy.zeros(dsize);
    for i in range(points.shape[0]):
        if points[i,0] > 0 and points[i,0] < dsize[0] and points[i,1] > 0 and points[i,1] < dsize[1] and points[i,2] > 0 and points[i,2] < dsize[2]:
               imgc[points[i,0], points[i,1], points[i,2]] = 1;
    self.plotOverlayLabel(image[:,:,0:dsize[2]], imgc, alpha = False);


if __name__ == "__main__":
    import numpy as np
    import iDISCO.Visualization.Plot as iplt
    l = np.array([[0,0,0,0,0], [0,1,1,0,0], [3,0,5,0,2], [5,0,0,0,0], [4,4,0,0,0]])
    x = np.random.rand(5,5);  
    
    iplt.plotOverlayLabel(x,l, alpha = False);    
    iplt.plotOverlayLabel(x,l, alpha = True);    
    

    
    
    
    
    

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

import iDISCO.IO.IO as io
import iDISCO.Analysis.Voxelization as vox


def plotTiling(image, tiling = "automatic", maxtiles = 20): 
    """Plot 3d image as tiles"""
    
    dim = image.ndim;
    
    if dim < 2 or dim > 4:
        raise StandardError('plotTiling: image dimension must be 2 to 4');    
    
    if dim == 2:
        image = image.reshape(image.shape + (1,));
        dim = 3;

    if image.ndim == 3: 
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
        imgpl = imgpl.transpose([1,0,2]);      
        
        if imgpl.shape[2] == 1:
            imgpl = imgpl.reshape((imgpl.shape[0], imgpl.shape[1]));       
        #a.imshow(imgpl, interpolation='none', cmap = cmap, vmin = imin, vmax = imax);
        a.imshow(imgpl, interpolation='none', cmap = cmap, vmin = imin, vmax = imax);
    
    fig.canvas.manager.window.activateWindow()
    fig.canvas.manager.window.raise_()
    
    return fig;


def overlayLabel(image, label, alpha = False, labelColorMap = 'jet'):
    """Overlay a gray scale image with colored labeled image"""
    
    lmax = label.max();
    
    if lmax <= 1:
        carray = numpy.array([[1,0,0,1]]);
    else:
        cm = mpl.cm.get_cmap(labelColorMap);
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
    return plotTiling(ov, tiling = tiling);



def overlayPoints(dataSource, pointSource, x = all, y = all, z = all, pointColor = [1,0,0]):
    """Overlay points on 3D Data"""
    data = io.readData(dataSource, x = x, y = y, z = z);
    points = io.readPoints(pointSource, x = x, y = y, z = z, shift = True);
    #print data.shape
    
    if not pointColor is None:
        dmax = data.max(); dmin = data.min();
        if dmin == dmax:
            dmax = dmin + 1;
        cimage = numpy.repeat( (data - dmin) / (dmax - dmin), 3);
        cimage = cimage.reshape(data.shape + (3,));    
    
        if data.ndim == 2:
            for p in points: # faster version using voxelize ?
                cimage[p[0], p[1], :] = pointColor;
        elif data.ndim == 3:
            for p in points: # faster version using voxelize ?
                cimage[p[0], p[1], p[2], :] = pointColor;
        else:
            raise RuntimeError('overlayPoints: data dimension %d not suported' % data.ndim);
    
    else:
        cimage = vox.voxelize(points, data.shape, voxelizationMethod = 'Pixel');
        cimage = cimage.astype(data.dtype) * data.max();
        data.shape = data.shape + (1,);
        cimage.shape =  cimage.shape + (1,);
        cimage = numpy.concatenate((data, cimage), axis  = 3);
    
    #print cimage.shape    
    return cimage;   


def plotOverlayPoints(dataSource, pointSource, pointColor = [1,0,0], x = all, y = all, z = all):
    """Plot points overlayed on gray scale 3d image"""
    cimg = self.overlayPoints(dataSource, pointSource, pointColor = pointColor, x = x, y = y, z = z);
    return self.plotTiling(cimg);
        



def test():
    """Test Plot module"""
    import numpy as np
    import iDISCO.Visualization.Plot as self
    reload(self)
    
    l = np.array([[0,0,0,0,0], [0,1,1,0,0], [3,0,5,0,2], [5,0,0,0,0], [4,4,0,0,0]])
    x = np.random.rand(5,5);  
    
    self.plotOverlayLabel(x,l, alpha = False);    
    self.plotOverlayLabel(x,l, alpha = True);    


    # 
    x = np.random.rand(50,20); 
    p = np.array([[10,15], [40,10]]);
    
    self.plotOverlayPoints(x, p)


if __name__ == "__main__":
    test();

    
    

    
    
    
    

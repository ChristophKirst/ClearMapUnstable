# -*- coding: utf-8 -*-

"""
Cell Size Detection Algorithms

Notes:
    -used for finding maxima and cell intensities

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import numpy

from scipy.ndimage.measurements import watershed_ift
#from skimage.morphology import watershed
from iDISCO.Analysis.Voxelization import voxelizePixel

#from skimage.measure import regionprops
import scipy.ndimage.measurements



from iDISCO.Utils.Timer import Timer
from iDISCO.Utils.ParameterTools import writeParameter;
from iDISCO.Visualization.Plot import plotOverlayLabel



import iDISCO.IO.IO as io


##############################################################################
# Basic Transforms
##############################################################################


def detectCellShape(img, peaks, cellShapeThreshold = 0, cellShapeFile = None, subStack = None, verbose = False, out = sys.stdout, **parameter):
    """Find cell shapes as labeled image"""  
    timer = Timer();cellShapeThreshold
    writeParameter(out = out, head = 'Cell Shape:', cellShapeThreshold = cellShapeThreshold);
    
    # extended maxima
    timer.reset();
    
    if cellShapeThreshold is None:
        imgmask = None;
    else:
        imgmask = img > cellShapeThreshold;
    
    imgpeaks = voxelizePixel(peaks, dataSize = img.shape, weights = numpy.arange(1, peaks.shape[0]+1));
    
    #imgws = watershed(img, imgpeaks, mask = imgmask);
    imgws = watershed_ift(img.astype('uint16'), imgpeaks);
    imgws[numpy.logical_not(imgmask)] = 0;
    
    if not cellShapeFile is None:
        if not subStack is None:
            ii = subStack["zSubStackCenterIndices"][0];
            ee = subStack["zSubStackCenterIndices"][1];
            si = subStack["zCenterIndices"][0];
        else:
            si = 0;
            ii = 0;
            ee = -1;
        
        io.writeData(cellShapeFile, imgws[:,:,ii:ee].astype('int32'), startIndex = si );
    
    
    if verbose:
        #plotTiling(img)
        plotOverlayLabel(img * 0.01, imgws, alpha = False);
        #plotOverlayLabel(img, imgmax.astype('int64'), alpha = True)     
    
    out.write(timer.elapsedTime(head = 'Cell Shape:') + '\n');
    return imgws


def findCellSize(imglabel, out = sys.stdout, **parameter):
    """Find cell size given cell shapes as labled image"""  
    timer = Timer();
    #writeParameter(out = out, head = 'Cell Size:');
        
    size = scipy.ndimage.measurements.sum(numpy.ones(imglabel.shape, dtype = bool), labels = imglabel, index = numpy.arange(1,imglabel.max() + 1));
        
    out.write(timer.elapsedTime(head = 'Cell Size:') + '\n');
    return size


def findCellIntensity(img, imglabel, cellIntensityMethod = 'Sum', out = sys.stdout, **parameter):
    """Find integrated cell intensity given cell shapes as labled image"""  
    timer = Timer();
    writeParameter(out = out, head = 'Cell Intensity:', cellIntensityMethod = cellIntensityMethod);
    
    if cellIntensityMethod == 'Sum':
        i = scipy.ndimage.measurements.sum(img, labels = imglabel, index = numpy.arange(1,imglabel.max() + 1));
    elif cellIntensityMethod == 'Mean':
        i = scipy.ndimage.measurements.mean(img, labels = imglabel, index = numpy.arange(1,imglabel.max() + 1));
    elif cellIntensityMethod == 'Max':
        i = scipy.ndimage.measurements.maximum(img, labels = imglabel, index = numpy.arange(1,imglabel.max() + 1));
    elif cellIntensityMethod == 'Min':
        i = scipy.ndimage.measurements.minimum(img, labels = imglabel, index = numpy.arange(1,imglabel.max() + 1));
    else:
        raise RuntimeError('cellIntensity: unkown method %s!' % cellIntensityMethod);
        
    out.write(timer.elapsedTime(head = 'Cell Intensity:') + '\n');
    
    return i
    
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:50:57 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import numpy
import pyximport;
#pyximport.install();
pyximport.install(setup_args={"include_dirs":numpy.get_include()}, reload_support=True)

import iDISCO.Analysis.VoxelizationCode as vox

import iDISCO.IO.IO as io
import math


def voxelize(points, dataSize = None, sink = None, voxelizationSize = (5,5,5), voxelizationMethod = 'Spherical', voxelizationWeights = None):
    
    if dataSize is None:
        dataSize = tuple(int(math.ceil(points[:,i].max())) for i in range(points.shape[1]));
    elif isinstance(dataSize, basestring):
        dataSize = io.dataSize(dataSize);
    
    points = io.readPoints(points);
        
    if voxelizationMethod == 'Spherical':
        if voxelizationWeights is None:
            data = vox.voxelizeSphere(points, dataSize[0], dataSize[1], dataSize[2], voxelizationSize[0], voxelizationSize[1], voxelizationSize[2]);
        else:
            data = vox.voxelizeSphereWithWeights(points, dataSize[0], dataSize[1], dataSize[2], voxelizationSize[0], voxelizationSize[1], voxelizationSize[2], voxelizationWeights);
           
    elif voxelizationMethod == 'Rectangular':
        if voxelizationWeights is None:
            data = vox.voxelizeRectangle(points, dataSize[0], dataSize[1], dataSize[2], voxelizationSize[0], voxelizationSize[1], voxelizationSize[2]);
        else:
            data = vox.voxelizeRectangleWithWeights(points, dataSize[0], dataSize[1], dataSize[2], voxelizationSize[0], voxelizationSize[1], voxelizationSize[2], voxelizationWeights);
    
    elif voxelizationMethod == "Pixel":
        data = self.voxelizePixel(points, dataSize, voxelizationWeights);
        
    else:
        raise RuntimeError('voxelize: mode: %s not supported!' % voxelizationMethod);
    
    return io.writeData(sink, data);


def voxelizePixel(points,  dataSize = None, weights = None):
    """Mark single pixel of each point in image"""
    
    if dataSize is None:
        dataSize = tuple(int(math.ceil(points[:,i].max())) for i in range(points.shape[1]));
    elif isinstance(dataSize, basestring):
        dataSize = io.dataSize(dataSize);
    
    if weights is None:
        vox = numpy.zeros(dataSize, dtype=numpy.int16);
        for i in range(points.shape[0]):
            if points[i,0] > 0 and points[i,0] < dataSize[0] and points[i,1] > 0 and points[i,1] < dataSize[1] and points[i,2] > 0 and points[i,2] < dataSize[2]:
                vox[points[i,0], points[i,1], points[i,2]] = 1;
    else:
        vox = numpy.zeros(dataSize, dtype=weights.dtype);
        for i in range(points.shape[0]):
            if points[i,0] > 0 and points[i,0] < dataSize[0] and points[i,1] > 0 and points[i,1] < dataSize[1] and points[i,2] > 0 and points[i,2] < dataSize[2]:
                vox[points[i,0], points[i,1], points[i,2]] = weights[i];
    
    return  vox;



def test():
    """Test vioxelization Module"""
    #import iDISCO.Analysis.Voxelization as self
    
    import iDISCO.Analysis.VoxelizationCode as vox
    import numpy
    
    
    points = numpy.random.rand(200,3) * 10;
    
    #use cython code
    vi = vox.voxelizeSphere(points, 20,20,20, 5,5,5);
    
    import iDISCO.Visualization.Plot as Plot

    Plot.plotTiling(vi)
    
    #use voxelize
    vi = self.voxelize(points, dataSize = (20,20,20), average = (5,5,5));
    
    Plot.plotTiling(vi)
    
    
    #weighted voxelization 
    
    points = numpy.random.rand(10,3) * 10;    
    weights = numpy.random.rand(10);
    
    #use voxelize
    vi = self.voxelize(points, dataSize = (20,20,20), average = (5,5,5));
    viw =  self.voxelize(points, dataSize = (20,20,20), average = (5,5,5), weights = weights);
    
    Plot.plotTiling(vi)
    Plot.plotTiling(viw)
    
   
if __name__ == "__main__":
    self.test();
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


def voxelize(points, imagesize, average = (5,5,5), mode = 'Spherical'):
    
    if mode == 'Spherical':
        return vox.voxelizeSphere(points, imagesize[0], imagesize[1], imagesize[2], average[0], average[1], average[2]);
    elif mode == 'Rectangular':
        return vox.voxelizeRectangle(points, imagesize[0], imagesize[1], imagesize[2], average[0], average[1], average[2]);   
    elif mode == "Pixel":
        return self.voxelizePixel(points, imagesize); 
    else:
        raise RuntimeError('voxelize: mode: %s not supported!' % mode);


def voxelizePixel(points, imagesize):
    vox = numpy.zeros(imagesize, dtype=numpy.int32);
    for i in range(points.shape[0]):
        if points[i,0] > 0 and points[i,0] < imagesize[0] and points[i,1] > 0 and points[i,1] < imagesize[1] and points[i,2] > 0 and points[i,2] < imagesize[2]:
               vox[points[i,0], points[i,1], points[i,2]] = 1;
    
    return  vox;



def test():
    points = numpy.random.rand(200,3) * 10;

    #use cython code
    vi = vox.voxelizeSphere(points, 20,20,20, 5,5,5);
    
    import iDISCO.Visualization.Plot as Plot

    Plot.plotTiling(vi)
    
    #use voxelize
    vi = self.voxelize(points, 20,20,20, 5,5,5);
    
    Plot.plotTiling(vi)
    
    
if __name__ == "__main__":
    self.test();
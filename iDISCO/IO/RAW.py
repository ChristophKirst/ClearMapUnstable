# -*- coding: utf-8 -*-
"""
Simple Interface to read RAW/MHD files as created by elastix

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import vtk
from vtk.util.numpy_support import vtk_to_numpy


def readData(filename, x = all, y = all, z = all, channel = 0, timepoint = 0, resolution = 0):
    """Read data from raw/mhd image"""
    
    imr = vtk.vtkMetaImageReader()
    imr.SetFileName(filename);
    imr.Update()

    im = imr.GetOutput()
    dims = im.GetDimensions()
    sc = im.GetPointData().GetScalars()
    img = vtk_to_numpy(sc)
    #print img.shape

    dims = (dims[2], dims[1], dims[0]);
    img = img.reshape(dims)
    img = img.transpose([1,2,0]);
    
    return img;


def test():    
    from iDISCO.Parameter import iDISCOPath
    import os
    
    """Test RAW module"""
    basedir = iDISCOPath();
    fn = os.path.join(basedir, 'Test/Synthetic/transformix/result.mhd')

    print "Loading rae image from: " + fn;
    
    img = self.readData(fn);  
    print "Image size: " + str(img.shape)
    
    import iDISCO.Visualization.Plot as Plot
    Plot.plotTiling(0.01 * img)
    
    #compare with original
    import iDISCO.IO.IO as io
    imgo = io.readData(os.path.join(basedir, 'Test/Synthetic/test_iDISCO_resample.tif'));
    Plot.plotTiling(0.01 * imgo)
     

if __name__ == "__main__":
    self.test();
    

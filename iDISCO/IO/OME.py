# -*- coding: utf-8 -*-
"""
Simple Interface to OME files

read data


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import 
import numpy

import sys
self = sys.modules[__name__];


def readData(filename, x = all, y = all, z = all):
    
    f = h5py.File(filename, "r");
    dataset = self.readDataSet(f, resolution=0);
    datasetsize = dataset.shape;
    
    if x == all:
        x = (0, datasetsize[0]);
    
    if y == all:
        y = (0, datasetsize[1]);
    
    if z == all:
        z = (0, datasetsize[2]);
        
    
    img = dataset[z[0]:z[1],y[0]:y[1], z[0]:z[1]];
    img = img.transpose((2,1,0)); # imaris stores files in reverse x,y,z ordering
    #img = dataset[1200:1400,1200:1400, zrange[0]:zrange[1]];   
    #img = dataset[0:50, 0:50, zrange[0]:zrange[1]];
    
    f.close();
    
    return img;

  



# -*- coding: utf-8 -*-
"""
Simple Interface to OME files

read data


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import numpy
from mahotas import imread

import sys
self = sys.modules[__name__];

import os
import re


def readData(filename, x = all, y = all, z = all):
    """Read data from individual image z slices"""
        
    #get path        
    (fpath, fname) = os.path.split(filename)
    fnames = os.listdir(fpath);
    
    searchRegex = re.compile(fname).search
    res= [ ( l, m.group(1) ) for l in fnames for m in (searchRegex(l),) if m]
    
    if res == []:
        raise RuntimeError('no files found in ' + fpath + ' match ' + fname + ' !');
        
    res.sort();
    nz = len(res);
    
    #determine zrange
    if z == all:
        #assume slices go from 0 to n-1
        z = (0, nz-1);
    else:
        nz = z[1] - z[0];

    fn = os.path.join(fpath, res[0]);
    img = imread(fn);
            
    if x== all:
        nx = img.shape[0];
        x = (0, nx);
    else:
        nx = x[1] - x[0];
            
    if y == all:
        ny = img.shape[0];
        y = (0, ny);
    else:
        ny = y[1] - y[0];        
    
    img = numpy.zeros((nx,ny,nz), dtype = img.dtype);
    
    for iz in range(z[0], z[1]):
        img[:,:,iz - z[0]] = imread(res[iz]);
    
    return img;

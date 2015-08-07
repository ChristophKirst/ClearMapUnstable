# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 19:25:15 2015

@author: mtllab
"""

import numpy
import iDISCO.IO.IO as io

import sys
self = sys.modules[__name__];

def labelPoints(points, labelimage):
    
    #points are (y,x,z)
    x = points[:,1];
    y = points[:,0];
    z = points[:,2];    
    nPoint = x.size;    
    
    pointLabels = numpy.ones(nPoint, 'int32');
    if isinstance(labelimage, basestring):
        labelImage = io.readData(labelimage);
                  
    dsize = labelImage.shape;
    for i in range(nPoint):
        if y[i] >= 0 and y[i] < dsize[0] and x[i] >= 0 and x[i] < dsize[1] and z[i] >= 0 and z[i] < dsize[2]:
             pointLabels[i] = labelImage[y[i], x[i], z[i]];
    
    return pointLabels;
    
    
def countPointsInRegions(points, labelimage):
    pointLabels = self.labelPoints(points, labelimage);
    
    cc = numpy.bincount(pointLabels);
    ii = numpy.nonzero(cc)[0];
    return numpy.vstack((ii,cc[ii])).T;
    
    

    
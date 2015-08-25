# -*- coding: utf-8 -*-
"""
Interface to write binary files of cell coordinates

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import numpy

import sys
self = sys.modules[__name__];

import iDISCO.IO.IO as io;

def writePoints(filename, points, **args):
    numpy.save(filename, points)
    return filename


def readPoints(filename, **args):
    points = numpy.load(filename);
    return io.pointsToRange(points, **args);


def test():    
    """Test NPY module"""
    
    import os, numpy
    import iDISCO.IO.NPY as self
    
    fn = os.path.split(self.__file__);
    fn = os.path.join(fn[0], '../Test/Data/NPY/points.npy');
    
    points = numpy.random.rand(5,3);
    self.writePoints(fn, points);  
    print "Wrote points to " + fn;
    print "Points:"
    print points
    
    points2 = self.readPoints(fn);
    print "Read points: "
    print points2
    
    print "Difference: " + str(numpy.abs(points-points2).max())
    

if __name__ == "__main__":
    
    self.test();
    

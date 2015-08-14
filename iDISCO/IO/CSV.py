# -*- coding: utf-8 -*-
"""
Simple Interface to write csv files of cell coordinates

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import numpy

import sys
self = sys.modules[__name__];

import iDISCO.IO.IO as io;

def writePoints(filename, points):
    numpy.savetxt(filename, points, delimiter=',', newline='\n', fmt='%.5e')

def readPoints(filename, **args):
    points = numpy.loadtxt(filename, delimiter=',');
    return io.pointsToRange(points, **args);
    

def test():    
    """Test CSV module"""
    
    import os
    
    fn = os.path.split(self.__file__);
    fn = os.path.join(fn[0], '../Test/ImageProcessing/points.txt');
    
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
    

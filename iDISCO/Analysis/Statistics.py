# -*- coding: utf-8 -*-
"""
Create some statistics to test significant changes
in voxelized and labeled data

Created on Mon Aug 10 16:04:25 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import numpy
from scipy import stats

import iDISCO.IO.IO as io


def readGroup(filenames, **args):
    """Turn a list of filenames into a numpy stack"""
    
    #check if stack already:
    if isinstance(filenames, numpy.ndarray):
        return filenames;
    
    #read the individual files
    group = [];
    for f in filenames:
        data = io.readData(f, **args);
        data = numpy.reshape(data, (1,) + data.shape);
        group.append(data);
    
    return numpy.vstack(group);


def tTest(group1, group2, signed = False):
    """t-Test on differences between the individual voxels in group1 and group2, group is a array of voxelizations"""
    
    g1 = self.readGroup(group1);  
    g2 = self.readGroup(group2);  
    
    tvals, pvals = stats.ttest_ind(g1, g2, axis = 0, equal_var = False);
    
    if signed:
        pvals = pvals * numpy.sign(tvals);
    
    return pvals;
    

def mean(group, **args):
    g = self.readGroup(group, **args);  
    return g.mean(axis = 0);


def std(group, **args):
    g = self.readGroup(group, **args);  
    return g.std(axis = 0);

   
def var(group, **args):
    g = self.readGroup(group, **args);  
    return g.var(axis = 0);    
    



    
def thresholdPoints(points, intensities, threshold = 0, row = 0):
    """Threshold points by intensities"""
    
    points, intensities = io.readPoints((points, intensities));   
    
    if intensities.ndim > 1:
        i = intensities[:,row];
    else:
        i = intensities;
        
    if not isinstance(threshold, tuple):
        threshold = (threshold, all);    
    
    iids = numpy.ones(i.shape, dtype = 'bool');
    if not threshold[0] is all:
        iids = numpy.logical_and(iids, i >= threshold[0]);
    if not threshold[1] is all:
        iids = numpy.logical_and(iids, i <= threshold[1]);
    
    return (points[iids, :], intensities[iids, :]);



def weightsFromPrecentiles(intensities, percentiles = [25,50,75,100]):
    perc = numpy.percentiles(intensities, percentiles);
    weights = numpy.zeros(intensities.shape);
    for p in perc:
        ii = intensities > p;
        weights[ii] = weights[ii] + 1;
    
    return weights;
        

def colorPValues(pvals, positive = [1,0,0], negative = [0,1,0]):
    ds = pvals.shape + (3,);
    pvc = numpy.zeros(ds);
    
    ids = pvals > 0;
    pvalsi = pvals[ids];
    for i in range(3):
        pvc[ids, i] = pvalsi * positive[i];
    
    ids = pvals < 0;
    pvalsi = pvals[ids];
    for i in range(3):
        pvc[ids, i] = -pvalsi * negative[i];
        
    return pvc;


def test():
    """Test the statistics array"""
    #import iDISCO.Analysis.Statistics as self
    # reload(self)
    #x = stats.norm.rvs(loc=5,scale=1,size=1500)
    #y = stats.norm.rvs(loc=-5,scale=1,size=1500)
    x = numpy.random.rand(4,3,20);
    y = numpy.random.rand(5, 3,20) + 0.5;
    
    print stats.ttest_ind(x,y, axis = 0, equal_var = False);
    
    pvals = self.tTest(x,y);
    
    print pvals
    
    
if __name__ == "__main__":
    self.test();
    
    

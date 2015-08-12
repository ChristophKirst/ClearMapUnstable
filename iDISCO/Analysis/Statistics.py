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



def tTest(group1, group2):
    """t-Test on differences between the individual voxels in group1 and group2, group is a array of voxelizations"""
    
    g1 = self.readGroup(group1);  
    g2 = self.readGroup(group2);  
    
    tvals, pvals = stats.ttest_ind(g1, g2, axis = 0, equal_var = False);
    
    return pvals;
    



def test():
    #import iDISCO.Analysis.Statistics as self
    # reload(self)
    #x = stats.norm.rvs(loc=5,scale=1,size=1500)
    #y = stats.norm.rvs(loc=-5,scale=1,size=1500)
    x = numpy.random.rand(15,3,20);
    y = numpy.random.rand(5, 3,20) + 0.5;
    
    print stats.ttest_ind(x,y, axis = 0, equal_var = False);
    
    pvals = self.tTest(x,y);
    
    print pvals
    
    
if __name__ == "__main__":
    self.test();
    
    

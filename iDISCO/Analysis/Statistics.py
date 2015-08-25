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


def tTest(group1, group2, signed = False, removeNaN = True, pcutoff = None):
    """t-Test on differences between the individual voxels in group1 and group2, group is a array of voxelizations"""
    
    g1 = self.readGroup(group1);  
    g2 = self.readGroup(group2);  
    
    tvals, pvals = stats.ttest_ind(g1, g2, axis = 0, equal_var = False);

    #remove nans
    if removeNaN: 
        pi = numpy.isnan(pvals);
        pvals[pi] = 1.0;
        tvals[pi] = 0;

    pvals = self.cutoffPValues(pvals, pcutoff = pcutoff);

    #return
    if signed:
        return pvals, numpy.sign(tvals);
    else:
        return pvals;
        
def cutoffPValues(pvals, pcutoff = 0.05):
    if pcutoff is None:
        return pvals;
    
    pvals2 = pvals.copy();
    pvals2[pvals2 > pcutoff]  = pcutoff;
    return pvals2;
    

def colorPValues(pvals, psign, positive = [1,0], negative = [0,1]):
    #invert pvals    
    pmax = pvals.max();
    pvalsinv = pvals.copy();
    pvalsinv = pmax - pvalsinv;
        
    d = len(positive);
    ds = pvals.shape + (d,);
    pvc = numpy.zeros(ds);
    
    #color
    ids = psign > 0;
    pvalsi = pvalsinv[ids];
    for i in range(d):
        pvc[ids, i] = pvalsi * positive[i];
    
    ids = psign < 0;
    pvalsi = pvalsinv[ids];
    for i in range(d):
        pvc[ids, i] = pvalsi * negative[i];
        
    return pvc;
    
    
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
            
    if not isinstance(threshold, tuple):
        threshold = (threshold, all);    
    
    if not isinstance(row, tuple):
        row = (row, row);
    
    
    if intensities.ndim > 1:
        i = intensities[:,row[0]];
    else:
        i = intensities;    
    
    iids = numpy.ones(i.shape, dtype = 'bool');
    if not threshold[0] is all:
        iids = numpy.logical_and(iids, i >= threshold[0]);
        
    if intensities.ndim > 1:
        i = intensities[:,row[1]];
    
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
        



def test():
    """Test the statistics array"""
    import iDISCO.Analysis.Statistics as self
    reload(self)
    import numpy
    #x = stats.norm.rvs(loc=5,scale=1,size=1500)
    #y = stats.norm.rvs(loc=-5,scale=1,size=1500)
    s = numpy.ones((5,4,20));
    s[:, 0:3, :] = - 1;
    
    x = numpy.random.rand(4,4,20);
    y = numpy.random.rand(5,4,20) + s;
    
    #print stats.ttest_ind(x,y, axis = 0, equal_var = False);
    
    pvals = self.tTest(x,y);
    
    print pvals
    
    pvals = self.tTest(x,y, signed = True);
    pvalscol = self.colorPValues(pvals)
    
    import iDISCO.Visualization.Plot as plt
    plt.plotTiling(1/pvalscol)
    
    
if __name__ == "__main__":
    self.test();
    
    

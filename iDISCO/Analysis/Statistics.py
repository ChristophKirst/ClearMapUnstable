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

import iDISCO.Analysis.Label as lbl


def readDataGroup(filenames, combine = True, **args):
    """Turn a list of filenames for data into a numpy stack"""
    
    #check if stack already:
    if isinstance(filenames, numpy.ndarray):
        return filenames;
    
    #read the individual files
    group = [];
    for f in filenames:
        data = io.readData(f, **args);
        data = numpy.reshape(data, (1,) + data.shape);
        group.append(data);
    
    if combine:
        return numpy.vstack(group);
    else:
        return group;


def readPointsGroup(filenames, **args):
    """Turn a list of filenames for points into a numpy stack"""
    
    #check if stack already:
    if isinstance(filenames, numpy.ndarray):
        return filenames;
    
    #read the individual files
    group = [];
    for f in filenames:
        data = io.readPoints(f, **args);
        #data = numpy.reshape(data, (1,) + data.shape);
        group.append(data);
    
    return group
    #return numpy.vstack(group);


def tTestVoxelization(group1, group2, signed = False, removeNaN = True, pcutoff = None):
    """t-Test on differences between the individual voxels in group1 and group2, group is a array of voxelizations"""
    
    g1 = self.readDataGroup(group1);  
    g2 = self.readDataGroup(group2);  
    
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
        

def countPointsGroupInRegions(pointGroup, withIds = True, labeledImage = lbl.DefaultLabeledImageFile):
     "Generates a table of counts for the various point datasets in pointGroup"""
     
     counts = [lbl.countPointsInRegions(p, labeledImage = labeledImage, sort = True, allIds = True) for p in pointGroup];
         
     if withIds:
         ids = counts[0][:,0];
         ids.shape = (1,) + ids.shape;
         counts = numpy.vstack((c[:,1] for c in counts)).T;
         return numpy.concatenate((ids.T,counts), axis = 1)
     else:
         return numpy.vstack((c[:,1] for c in counts)).T;
         

def tTestPointsInRegions(pointGroup1, pointGroup2, labeledImage = lbl.DefaultLabeledImageFile, withIds = True, signed = False, removeNaN = True, pcutoff = None):
    """t-Test on differences in counts of points in labeled regions"""
    
    p1 = countPointsGroupInRegions(pointGroup1, labeledImage = labeledImage);
    ids = p1[:,0];
    
    p1  = p1[:,1::];
     
    p2 = countPointsGroupInRegions(pointGroup1,  labeledImage = labeledImage, withIds = False);   
    
    tvals, pvals = stats.ttest_ind(p1, p2, axis = 1, equal_var = False);

    #remove nans
    if removeNaN: 
        pi = numpy.isnan(pvals);
        pvals[pi] = 1.0;
        tvals[pi] = 0;

    pvals = self.cutoffPValues(pvals, pcutoff = pcutoff);

    if withIds:
        pvals.shape = (1,) + pvals.shape;
        ids.shape = (1,) + ids.shape;
        pvals = numpy.concatenate((ids.T, pvals.T), axis = 1);
        
    #return
    if signed:
        return pvals, numpy.sign(tvals);
    else:
        return pvals;

    


def test():
    """Test the statistics array"""
    import iDISCO.Analysis.Statistics as self
    reload(self)
    import numpy, os
    #x = stats.norm.rvs(loc=5,scale=1,size=1500)
    #y = stats.norm.rvs(loc=-5,scale=1,size=1500)
    s = numpy.ones((5,4,20));
    s[:, 0:3, :] = - 1;
    
    x = numpy.random.rand(4,4,20);
    y = numpy.random.rand(5,4,20) + s;
    
    # print stats.ttest_ind(x,y, axis = 0, equal_var = False);
    pvals, psign = self.tTestVoxelization(x,y, signed = True);
    print pvals
    
    pvalscol = self.colorPValues(pvals, psign, positive = [255,0,0], negative = [0,255,0])
    
    import iDISCO.Visualization.Plot as plt
    plt.plotTiling(pvalscol)
    
    # test points
    import iDISCO.Settings as settings
    pf = os.path.join(settings.IDISCOPath, 'Test/Synthetic/cells_transformed_to_reference.csv');
    
    pg = (pf,pf);
    
    pc = self.countPointsGroupInRegions(pg);
        
    pvals, tvals = self.tTestPointsInRegions(pg, pg, signed = True);
    
    
if __name__ == "__main__":
    self.test();
    
    

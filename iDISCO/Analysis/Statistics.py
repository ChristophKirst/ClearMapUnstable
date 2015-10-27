# -*- coding: utf-8 -*-
"""
Create some statistics to test significant changes
in voxelized and labeled data

Created on Mon Aug 10 16:04:25 2015

TODO; cleanup / make generic

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import numpy
from scipy import stats

import iDISCO.IO.IO as io
import iDISCO.Analysis.Label as lbl
import iDISCO.Analysis.StatisticsUtils as stats2

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
    
    tvals, pvals = stats.ttest_ind(g1, g2, axis = 0, equal_var = True);

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
        

# needs clean up
def countPointsGroupInRegions(pointGroup, labeledImage = lbl.DefaultLabeledImageFile, intensityGroup = None, intensityRow = 0, returnIds = True, returnCounts = False):
     """Generates a table of counts for the various point datasets in pointGroup"""
 
     if intensityGroup is None: 
         counts = [lbl.countPointsInRegions(pointGroup[i], labeledImage = labeledImage, sort = True, allIds = True, returnIds = False, returnCounts = returnCounts, intensities = None) for i in range(len(pointGroup))];
     else:
         counts = [lbl.countPointsInRegions(pointGroup[i], labeledImage = labeledImage, sort = True, allIds = True, returnIds = False, returnCounts = returnCounts,
                                            intensities = intensityGroup[i], intensityRow = intensityRow) for i in range(len(pointGroup))];
     
     if returnCounts and not intensityGroup is None:
         countsi = (c[1] for c in counts);
         counts  = (c[0] for c in counts);
     else:
         countsi = None;
         
     counts = numpy.vstack((c for c in counts)).T;
     if not countsi is None:
         countsi =  numpy.vstack((c for c in countsi)).T;
     
     if returnIds:
         ids = numpy.sort(lbl.Label.ids);
         #ids.shape = (1,) + ids.shape;
         
         #return numpy.concatenate((ids.T,counts), axis = 1
         if countsi is None:
             return ids, counts
         else:
             return ids, counts, countsi      
     else:
         if countsi is None:
             return counts
         else:
             return counts, countsi
         


# needs clean up
def tTestPointsInRegions(pointCounts1, pointCounts2, labeledImage = lbl.DefaultLabeledImageFile, signed = False, removeNaN = True, pcutoff = None, equal_var = False):
    """t-Test on differences in counts of points in labeled regions"""
    
    #ids, p1 = countPointsGroupInRegions(pointGroup1, labeledImage = labeledImage, withIds = True);
    #p2 = countPointsGroupInRegions(pointGroup2,  labeledImage = labeledImage, withIds = False);   
    
    tvals, pvals = stats.ttest_ind(pointCounts1, pointCounts2, axis = 1, equal_var = equal_var);
    
    #remove nans
    if removeNaN: 
        pi = numpy.isnan(pvals);
        pvals[pi] = 1.0;
        tvals[pi] = 0;

    pvals = self.cutoffPValues(pvals, pcutoff = pcutoff);
    
    #pvals.shape = (1,) + pvals.shape;
    #ids.shape = (1,) + ids.shape;
    #pvals = numpy.concatenate((ids.T, pvals.T), axis = 1);
    
    if signed:
        return pvals, numpy.sign(tvals);
    else:
        return pvals;

    


def testCompletedCumulatives(data, method = 'AndersonDarling', offset = None):
    """Test if data sets have the same number / intensity distribution by adding max intensity counts to the smaller sized data sets and performing a distribution comparison test"""
    
    #idea: fill up data points to the same numbers at the high intensity values and use KS test
    #cf. work in progress on thoouroghly testing the differences in histograms
    
    #fill up the low count data
    n = numpy.array([x.size for x in data]);
    nm = n.max();
    m = numpy.array([x.max() for x in data]);
    mm = m.max();
    k = n.size;
    #print nm, mm, k
    
    if offset is None:
        #assume data starts at 0 !
        offset = mm / nm; #ideall for all statistics this should be mm + eps to have as little influence as possible.
    

    datac = [x.copy() for x in data];
    for i in range(m.size):
        if n[i] < nm:
            datac[i] = numpy.concatenate((datac[i], numpy.ones(nm-n[i], dtype = datac[i].dtype) * (mm + offset)));
         
    #test by plotting
    import matplotlib.pyplot as plt;
    for i in range(m.size):
        datac[i].sort();
        plt.step(datac[i], numpy.arange(datac[i].size));
    
    #perfomr the tests
    if method == 'KolmogorovSmirnov' or method == 'KS':
        if k == 2:
            (s, p) = stats.ks_2samp(datac[0], datac[1]);
        else:
            raise RuntimeError('KolmogorovSmirnov only for 2 samples not %d' % k);
        
    elif method == 'CramervonMises' or method == 'CM':
        if k == 2:
            (s,p) = stats2.testCramerVonMises2Sample(datac[0], datac[1]);
        else:
            raise RuntimeError('CramervonMises only for 2 samples not %d' % k);
      
    elif method == 'AndersonDarling' or method == 'AD':
        (s,a,p) = stats.anderson_ksamp(datac);

    return (p,s);



def testCompletedCumulativesInSpheres(points1, points2, dataSize = lbl.DefaultLabeledImageFile, radius = 100, method = 'AndresonDarling'):
    """Performs completed cumulative distribution tests for each pixel using points in a ball centered at that cooridnates, returns 4 arrays p value, statistic value, number in each group"""
    
    #TODO: sinple implementation -> slow -> speed up
    dataSize = io.dataSize(dataSize);
    if len(dataSize) != 3:
        raise RuntimeError('dataSize expected to be 3d');
    
    # distances^2 to origin
    x1= points1[:,0]; y1 = points1[:,1]; z1 = points1[:,2]; i1 = points1[:,3];
    d1 = x1 * x1 + y1 * y1 + z1 * z1;
    
    x2 = points2[:,0]; y2 = points2[:,1]; z2 = points2[:,2]; i2 = points2[:,3];
    d2 = x2 * x2 + y2 * y2 + z2 * z2;
        
    r2 = radius * radius; # TODO: inhomogenous in 3d !
    
    p = numpy.zeros(dataSize);
    s = numpy.zeros(dataSize);
    n1 = numpy.zeros(dataSize, dtype = 'int');
    n2 = numpy.zeros(dataSize, dtype = 'int');
    
    for x in range(dataSize[0]):
        for y in range(dataSize[1]):
            for z in range(dataSize[2]):
                d11 = d1 - 2 (x * x1 + y * y1 + z * z1) + (x*x + y*y + z*z);
                d22 = d2 - 2 (x * x2 + y * y2 + z * z2) + (x*x + y*y + z*z);
                
                ii1 = d11 < r2;
                ii2 = d22 < r2;
                
                (p[x,y,z], s[x,y,z]) = self.testCompletedCumulatives((i1[ii1], i2[ii2]), method = method);
                n1[x,y,z] = ii1.sum();
                n2[x,y,z] = ii2.sum();
    
    return (p,s,n1,n2);
        


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
    
    

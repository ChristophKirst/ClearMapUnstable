# -*- coding: utf-8 -*-
"""
Process al large stack in parallel or sequentially

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import math
import numpy

from multiprocessing import Pool

import iDISCO.IO.IO as io

from iDISCO.ImageProcessing.SpotDetection import detectCells;

from iDISCO.Utils.ProcessWriter import ProcessWriter;
from iDISCO.Utils.Timer import Timer;

#define the subroutine for the processing

def processSubStack(dsr):
    """Helper to process stack in parallel"""
    
    fn = dsr[0];
    sf = dsr[1];
    pp = dsr[2];
    zr = dsr[3];
    xr = dsr[4];
    yr = dsr[5];
    ii = dsr[6];
    n  = dsr[7];
    
    pw = ProcessWriter(ii);
    timer = Timer();
    
    pw.write("processing chunk " + str(ii) + "/" + str(n));
    pw.write("file          = " + fn);
    pw.write("segmentation  = " + str(sf));
    pw.write("ranges: x,y,z = " + str(xr) +  "," + str(yr) + "," + str(zr)); 
    
    timer.reset();    
    img = io.readData(fn, z = zr, x = xr, y = yr); 
    print "hkjhkjhkjhkjh"
    print img.shape
    pw.write(timer.elapsedTime(head = 'Reading data of size ' + str(img.shape)));
    
    timer.reset();
    seg = sf(img, out = pw, **pp);
    pw.write(timer.elapsedTime(head = 'Segmenting data of size ' + str(img.shape)));  
    
    return seg;


def calculateChunkSize(size, processes = 2, chunkSizeMax = 100, chunkSizeMin = 30, chunkOverlap = 15,  chunkOptimization = True, chunkOptimizationSize = all, verbose = True):
    """Calculates the chunksize for parallel processing"""
    
    pre = "ChunkSize: ";
    
    #calcualte chunk sizes
    chunksize = chunkSizeMax;
    nchunks = int(math.ceil((size - chunksize) / (1. * (chunksize - chunkOverlap)) + 1)); 
    if nchunks <= 0:
        nchunks = 1;   
    chunksize = (size + (nchunks-1) * chunkOverlap) / nchunks;
    
    if verbose:
        print pre + "Estimated chunk size " + str(chunksize) + " in " + str(nchunks) + " chunks!";
    
    if nchunks == 1:
        return 1, [(0, chunksize)], [0, chunksize]
        
    #optimize number of chunks wrt to number of processors
    if chunkOptimization:
        np = nchunks % processes;
        if np != 0:
            if chunkOptimizationSize == all:
                if np < processes / 2.0:
                    chunkOptimizationSize = True;
                else:
                    chunkOptimizationSize = False;
                    
            if verbose:
                print pre + "Optimizing chunk size to fit number of processes!"
                
            if not chunkOptimizationSize:
                #try to deccrease chunksize / increase chunk number to fit distribution on processors
                nchunks = nchunks - np + processes;
                chunksize = (size + (nchunks-1) * chunkOverlap) / nchunks;
                
                if verbose:
                    print pre + "Optimized chunk size decreased to " + str(chunksize) + " in " + str(nchunks) + " chunks!";
                    
            else:
                if nchunks != np:
                    #try to decrease chunksize to have 
                    nchunks = nchunks - np;
                    chunksize = (size + (nchunks-1) * chunkOverlap) / nchunks;
                                  
                    if verbose:
                        print pre + "Optimized chunk size increased to " + str(chunksize) + " in " + str(nchunks) + " chunks!";
                
                else:
                    if verbose:
                        print pre + "Optimized chunk size unchanged " + str(chunksize) + " in " + str(nchunks) + " chunks!";
        
        else:
            if verbose:
                print pre + "Optimized chunk size unchanged " + str(chunksize) + " in " + str(nchunks) + " chunks!";
    
    
    #increase overlap if chunks to small
    chunkSizeMin = min(chunkSizeMin, chunkOverlap);
    if chunksize < chunkSizeMin:
        if verbose: 
            print pre + "Warning: optimal chunk size " + str(chunksize) + " smaller than minimum chunk size " + str(chunkSizeMin) + "!"; 
        chunksize = chunkSizeMin;
        chunkOverlap = math.ceil(chunksize - (size - chunksize) / (nchunks -1));
        
        if verbose:        
            print pre + "Warning: setting chunk overlap to " + str(chunkOverlap) + "!";
           
    #calucalte actual chunk sizes
    chunksizerest = chunksize;
    chunksize = int(math.floor(chunksize));
    chunksizerest = chunksizerest - chunksize;
    
    zranges = [(0, chunksize)];
    zcenters = [0];
    n = 1;
    csr = chunksizerest;
    zhi = chunksize;
    
    while (n < nchunks):
        n += 1;
        
        zhiold = zhi;
        zlo = zhi - chunkOverlap;
        zhi = zlo + chunksize;
        
        csr += chunksizerest;
        if csr >= 1:
            csr = csr - 1;
            zhi += 1;
        
        if n == nchunks:        
            zhi = size;
        
        zranges.append((int(zlo), int(zhi)));
        zcenters.append((zhiold - zlo) / 2. + zlo); 
        
    zcenters.append(size);
          
    
    if verbose:    
        print pre + "final chunks : " + str(zranges);
        print pre + "final centers: " + str(zcenters);
    
    return nchunks, zranges, zcenters;



def joinPoints(pointlist, zranges, zcenters, intensities = None, dataSize = None, x = all, y = all, z = all, shiftPoints = True, **args):
    """Joins a list of points obtained from processing a stack in chunks"""
    
    nchunks = len(pointlist);
    results = [];
    resultsi = [];
    for i in range(nchunks):
        cts = pointlist[i];
        if not intensities is None:
            cti = intensities[i];
        
        if cts.size > 0:
            cts[:,2] += zranges[i][0];
            iid = numpy.logical_and(cts[:,2] < zcenters[i+1], cts[:,2] >= zcenters[i]);
            cts = cts[iid,:];
            results.append(cts);
            if not intensities is None:
                cti = cti[iid];
                resultsi.append(cti);
            
    if results == []:
        if not intensities is None:
            return (numpy.zeros((0,3)), numpy.zeros((0)));
        else:
            return numpy.zeros((0,3))
    else:
        points = numpy.concatenate(results);
        
        if shiftPoints:
            points = points + io.pointShiftFromRange(dataSize, x = x, y = y, z = z);
            
        if intensities is None:
            return points;
        else:
            return (points, numpy.concatenate(resultsi));


def parallelProcessStack(source, x = all, y = all, z = all, sink = None,
                         processes = 2, chunkSizeMax = 100, chunkSizeMin = 30, chunkOverlap = 15,
                         chunkOptimization = True, chunkOptimizationSize = all, 
                         function = detectCells, join = self.joinPoints, **parameter):
    """Parallel process a stack using"""
    
    #determine z ranges
    fs = io.dataSize(source);
    zs = fs[2];
    zr = io.toDataRange(zs, r = z);
    nz = zr[1] - zr[0];
    nchunks, zranges, zcenters = self.calculateChunkSize(nz, processes = processes, chunkSizeMax = chunkSizeMax, chunkSizeMin = chunkSizeMin, chunkOverlap = chunkOverlap, chunkOptimization = chunkOptimization, chunkOptimizationSize = chunkOptimizationSize);
        
    #adjust for the zrange
    zcenters = [xc + zr[0] for xc in zcenters];
    zranges = [(zc[0] + zr[0], zc[1] + zr[0]) for zc in zranges];
    
    argdata = [];
    for i in range(nchunks):
        argdata.append((source, function, parameter, zranges[i], x, y, i, nchunks));
        
    #print argdata
    
    # process in parallel
    pool = Pool(processes = processes);
    #print processes
    
    resultsseg = pool.map(processSubStack, argdata);
    
    #print '=========== results';
    #print resultsseg;
    #print nchunks
        
    #join the results: extension: make this a join function and pass as argument 
    results = join([resultsseg[i][0] for i in range(nchunks)], zranges, zcenters, [resultsseg[i][1] for i in range(nchunks)], dataSize = fs, x = x, y = y, z = z, **parameter);
    
    #write / or return 
    return io.writePoints(sink, results);


def sequentiallyProcessStack(source, x = all, y = all, z = all, sink = None,
                             chunkSizeMax = 100, chunkSizeMin = 30, chunkOverlap = 15,
                             function = detectCells, join = self.joinPoints, **parameter):
    """Sequential image processing on a stack"""
    
    #determine z ranges
    fs = io.dataSize(source);
    zs = fs[2];
    zr = io.toDataRange(zs, r = z);
    nz = zr[1] - zr[0];
    nchunks, zranges, zcenters = self.calculateChunkSize(nz, processes = 1, chunkSizeMax = chunkSizeMax, chunkSizeMin = chunkSizeMin, chunkOverlap = chunkOverlap, chunkOptimization = False, chunkOptimizationSize = False);
        
    #adjust for the zrange
    zcenters = [xc + zr[0] for xc in zcenters];
    zranges = [(zc[0] + zr[0], zc[1] + zr[0]) for zc in zranges];
    
    argdata = [];
    for i in range(nchunks):
        argdata.append((source, function, parameter, zranges[i], x, y, i, nchunks));
        
    #run sequentially
    resultsseg = [];
    for i in range(nchunks):
        resultsseg.append(processSubStack(argdata[i]));
            
    #join the results  
    #join the results: extension: make this a join function and pass as argument 
    results = join([resultsseg[i][0] for i in range(nchunks)], zranges, zcenters, [resultsseg[i][1] for i in range(nchunks)], dataSize = fs, x = x, y = y, z = z, **parameter);

    #write / or return 
    return io.writePoints(sink, results);


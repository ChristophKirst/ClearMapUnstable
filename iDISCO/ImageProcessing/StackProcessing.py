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

from iDISCO.Parameter import ImageProcessingParameter
from iDISCO.ImageProcessing.SpotDetection import detectCells;


from iDISCO.Utils.ProcessWriter import ProcessWriter;
from iDISCO.Utils.Timer import Timer;

#define the subroutine for the processing

def processSubStack(dsr):
    
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
    pw.write(timer.elapsedTime(head = 'Reading data of size ' + str(img.shape)));
    
    timer.reset();
    seg = sf(img, out = pw, parameter = pp);
    pw.write(timer.elapsedTime(head = 'Segmenting data of size ' + str(img.shape)));  
    
    return seg;


def calculateChunkSize(size, chunksizemax = 100, chunksizemin = 30, chunkoverlap = 15, processes = 2, optimizechunks = True, optimizechunksizeincrease = all, verbose = True):
    """Calculates the chunksize for parallel processing"""
    
    pre = "ChunkSize: ";
    
    #calcualte chunk sizes
    chunksize = chunksizemax;
    nchunks = int(math.ceil((size - chunksize) / (1. * (chunksize - chunkoverlap)) + 1)); 
    if nchunks <= 0:
        nchunks = 1;   
    chunksize = (size + (nchunks-1) * chunkoverlap) / nchunks;
    
    if verbose:
        print pre + "Estimated chunk size " + str(chunksize) + " in " + str(nchunks) + " chunks!";
    
    if nchunks == 1:
        return 1, [(0, chunksize)], [0, chunksize]
        
    #optimize number of chunks wrt to number of processors
    if optimizechunks:
        np = nchunks % processes;
        if np != 0:
            if optimizechunksizeincrease == all:
                if np < processes / 2.0:
                    optimizechunksizeincrease = True;
                else:
                    optimizechunksizeincrease = False;
                    
            if verbose:
                print pre + "Optimizing chunk size to fit number of processes!"
                
            if not optimizechunksizeincrease:
                #try to deccrease chunksize / increase chunk number to fit distribution on processors
                nchunks = nchunks - np + processes;
                chunksize = (size + (nchunks-1) * chunkoverlap) / nchunks;
                
                if verbose:
                    print pre + "Optimized chunk size decreased to " + str(chunksize) + " in " + str(nchunks) + " chunks!";
                    
            else:
                if nchunks != np:
                    #try to decrease chunksize to have 
                    nchunks = nchunks - np;
                    chunksize = (size + (nchunks-1) * chunkoverlap) / nchunks;
                                  
                    if verbose:
                        print pre + "Optimized chunk size increased to " + str(chunksize) + " in " + str(nchunks) + " chunks!";
                
                else:
                    if verbose:
                        print pre + "Optimized chunk size unchanged " + str(chunksize) + " in " + str(nchunks) + " chunks!";
        
        else:
            if verbose:
                print pre + "Optimized chunk size unchanged " + str(chunksize) + " in " + str(nchunks) + " chunks!";
    
    
    #increase overlap if chunks to small
    chunksizemin = min(chunksizemin, chunkoverlap);
    if chunksize < chunksizemin:
        if verbose: 
            print pre + "Warning: optimal chunk size " + str(chunksize) + " smaller than minimum chunk size " + str(chunksizemin) + "!"; 
        chunksize = chunksizemin;
        chunkoverlap = math.ceil(chunksize - (size - chunksize) / (nchunks -1));
        
        if verbose:        
            print pre + "Warning: setting chunk overlap to " + str(chunkoverlap) + "!";
           
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
        zlo = zhi - chunkoverlap;
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



def parallelProcessStack(filename, x = all, y = all, z = all, 
                         processes = 2, chunksizemax = 100, chunksizemin = 30, chunkoverlap = 15,
                         optimizechunks = True, optimizechunksizeincrease = all, 
                         segmentation = detectCells, parameter = ImageProcessingParameter()):
    """Parallel segmetation on a stack distributes cell detection on parallel processes"""
    
    #determine z ranges
    zs = io.dataZSize(filename);
    zr = io.toDataRange(zs, r = z);
    nz = zr[1] - zr[0];
    nchunks, zranges, zcenters = self.calculateChunkSize(nz, chunksizemax = chunksizemax, chunksizemin = chunksizemin, chunkoverlap = chunkoverlap, processes = processes, optimizechunks = optimizechunks, optimizechunksizeincrease = optimizechunksizeincrease);
        
    #adjust for the zrange
    zcenters = [xc + zr[0] for xc in zcenters];
    zranges = [(zc[0] + zr[0], zc[1] + zr[0]) for zc in zranges];
    
    argdata = [];
    for i in range(nchunks):
        argdata.append((filename, segmentation, parameter, zranges[i], x, y, i, nchunks));
        
    #print argdata
    
    # process in parallel
    pool = Pool(processes = processes);
    #print processes
    
    resultsseg = pool.map(processSubStack, argdata);
    
    #print '=========== results';
    #print resultsseg;
    #print nchunks
        
    #join the results  
    results = [];
    resultsi = [];
    for i in range(nchunks):
        cts = resultsseg[i][0];
        cti = resultsseg[i][1];
        
        if cts.size > 0:
            cts[:,2] += zranges[i][0];
            iid = numpy.logical_and(cts[:,2] < zcenters[i+1], cts[:,2] >= zcenters[i]);
            cts = cts[iid,:];
            cti = cti[iid];
            results.append(cts);
            resultsi.append(cti);
            
    if results == []:
        return numpy.zeros((0,3)), numpy.zeros((0));
    else:
        return (numpy.concatenate(results), numpy.concatenate(resultsi));
        
        
        

def sequentiallyProcessStack(filename, x = all, y = all, z = all, 
                             chunksizemax = 100, chunksizemin = 30, chunkoverlap = 15,
                             segmentation = detectCells, parameter = ImageProcessingParameter()):
    """Sequential image processing on a stack"""
    
    #determine z ranges
    zr = io.readZRange(filename, z = z);
    nz = zr[1] - zr[0];
    nchunks, zranges, zcenters = self.calculateChunkSize(nz, chunksizemax = chunksizemax, chunksizemin = chunksizemin, chunkoverlap = chunkoverlap, processes = 1, optimizechunks = False, optimizechunksizeincrease = False);
        
    #adjust for the zrange
    zcenters = [xc + zr[0] for xc in zcenters];
    zranges = [(zc[0] + zr[0], zc[1] + zr[0]) for zc in zranges];
    
    argdata = [];
    for i in range(nchunks):
        argdata.append((filename, segmentation, parameter, zranges[i], x, y, i, nchunks));
        
    #run sequentially
    resultsseg = [];
    for i in range(nchunks):
        resultsseg.append(processSubStack(argdata[i]));
            
    #join the results  
    results = [];
    resultsi = [];
    for i in range(nchunks):
        cts = resultsseg[i][0];
        cti = resultsseg[i][1];
        
        if cts.size > 0:
            cts[:,2] += zranges[i][0];
            iid = numpy.logical_and(cts[:,2] < zcenters[i+1], cts[:,2] >= zcenters[i]);
            cts = cts[iid,:];
            cti = cti[iid];
            results.append(cts);
            resultsi.append(cti);
            
    if results == []:
        return numpy.zeros((0,3)), numpy.zeros((0));
    else:
        return (numpy.concatenate(results), numpy.concatenate(resultsi));
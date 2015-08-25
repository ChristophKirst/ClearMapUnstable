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


from iDISCO.Utils.ParameterTools import writeParameter
from iDISCO.Utils.ProcessWriter import ProcessWriter;
from iDISCO.Utils.Timer import Timer;

   
def printSubStackInfo(subStack, out = sys.stdout):
    writeParameter(head = "Sub Stack: ", out = out, **subStack);
    out.write('\n');


#define the subroutine for the processing
def processSubStack(dsr):
    """Helper to process stack in parallel"""

    sf  = dsr[0];
    pp  = dsr[1];
    sub = dsr[2];
    
    pw = ProcessWriter(sub["stackId"]);
    timer = Timer();
    
    pw.write("processing chunk " + str(sub["stackId"]) + "/" + str(sub["nStacks"]));
    pw.write("file          = " + sub["source"]);
    pw.write("segmentation  = " + str(sf));
    pw.write("ranges: x,y,z = " + str(sub["x"]) +  "," + str(sub["y"]) + "," + str(sub["z"])); 
    
    timer.reset();    
    img = io.readData(sub["source"], x = sub["x"], y = sub["y"], z = sub["z"]); 
    pw.write(timer.elapsedTime(head = 'Reading data of size ' + str(img.shape)));
    
    timer.reset();
    seg = sf(img, subStack = sub, out = pw, **pp);    
    
    pw.write(timer.elapsedTime(head = 'Segmenting data of size ' + str(img.shape)));
    
    return seg;


def joinPoints(results, subStacks = None, shiftPoints = True, **args):
    """Joins a list of points obtained from processing a stack in chunks"""
    
    nchunks = len(results);
    pointlist = [results[i][0] for i in range(nchunks)];
    intensities = [results[i][1] for i in range(nchunks)]; 
    
    results = [];
    resultsi = [];
    for i in range(nchunks):
        cts = pointlist[i];
        if not intensities is None:
            cti = intensities[i];
        
        if cts.size > 0:
            cts[:,2] += subStacks[i]["z"][0];
            iid = numpy.logical_and(cts[:,2] < subStacks[i]["zCenters"][1], cts[:,2] >=  subStacks[i]["zCenters"][0]);
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
            points = points + io.pointShiftFromRange(io.dataSize(subStacks[0]["source"]), x = subStacks[0]["x"], y = subStacks[0]["y"], z = 0);
        else:
            points = points - io.pointShiftFromRange(io.dataSize(subStacks[0]["source"]), x = 0, y = 0, z = subStacks[0]["z"]); #absolute offset is added initially via zranges !
            
        if intensities is None:
            return points;
        else:
            return (points, numpy.concatenate(resultsi));



def calculateChunkSize(size, processes = 2, chunkSizeMax = 100, chunkSizeMin = 30, chunkOverlap = 15,  chunkOptimization = True, chunkOptimizationSize = all, verbose = True):
    """Calculates the chunksize and other info for parallel processing, returns list of SubStack objects"""
    
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


def calculateSubStacks(source, z = all, x = all, y = all, **args):
    
    #determine z ranges
    fs = io.dataSize(source);
    zs = fs[2];
    zr = io.toDataRange(zs, r = z);
    nz = zr[1] - zr[0];
    
    #calculate optimal chunk sizes
    nchunks, zranges, zcenters = self.calculateChunkSize(nz, **args);
    
    #adjust for the zrange
    zcenters = [c + zr[0] for c in zcenters];
    zranges = [(zc[0] + zr[0], zc[1] + zr[0]) for zc in zranges];
    
    #create substacks
    subStacks = [];
    indexlo = zr[0];
    
    for i in range(nchunks):
        
        indexhi = int(round(zcenters[i+1]));
        if indexhi > zr[1] or i == nchunks - 1:
            indexhi = zr[1];
        
        zs = zranges[i][1] - zranges[i][0];
        
        subStacks.append({"stackId" : i, "nStacks" : nchunks, 
                          "source" : source, "x" : x, "y" : y, "z" : zranges[i], 
                          "zCenters" : (zcenters[i], zcenters[i+1]),
                          "zCenterIndices" : (indexlo, indexhi),
                          "zSubStackCenterIndices" : (indexlo - zranges[i][0], zs - (zranges[i][1] - indexhi))});
        
        indexlo = indexhi; # + 1;
    
    return subStacks;
        


def parallelProcessStack(source, x = all, y = all, z = all, sink = None,
                         processes = 2, chunkSizeMax = 100, chunkSizeMin = 30, chunkOverlap = 15,
                         chunkOptimization = True, chunkOptimizationSize = all, 
                         function = detectCells, join = self.joinPoints, **parameter):
    """Parallel process a stack"""
    
    subStacks = calculateSubStacks(source, x = x, y = y, z = z, 
                                   processes = processes, chunkSizeMax = chunkSizeMax, chunkSizeMin = chunkSizeMin, chunkOverlap = chunkOverlap,
                                   chunkOptimization = chunkOptimization, chunkOptimizationSize = chunkOptimizationSize, verbose = True);
                                   
    nSubStacks = len(subStacks);
    print "Number of SubStacks: %d" % nSubStacks;
                                       
    #for i in range(nSubStacks):
    #    self.printSubStackInfo(subStacks[i]);
    
    argdata = [];
    for i in range(nSubStacks):
        argdata.append((function, parameter, subStacks[i]));    
    #print argdata
    
    # process in parallel
    pool = Pool(processes = processes);    
    results = pool.map(processSubStack, argdata);
    
    #print '=========== results';
    #print results;
        
    #join the results
    results = join(results, subStacks = subStacks, **parameter);
    
    #write / or return 
    return io.writePoints(sink, results);


def sequentiallyProcessStack(source, x = all, y = all, z = all, sink = None,
                             chunkSizeMax = 100, chunkSizeMin = 30, chunkOverlap = 15,
                             function = detectCells, join = self.joinPoints, **parameter):
    """Sequential image processing on a stack"""
    
    #determine z ranges      
    subStacks = calculateSubStacks(source, x = x, y = y, z = z, 
                                   processes = 1, chunkSizeMax = chunkSizeMax, chunkSizeMin = chunkSizeMin, chunkOverlap = chunkOverlap,  
                                   chunkOptimization = False, verbose = True);
    
    nSubStacks = len(subStacks);
    argdata = [];
    for i in range(nSubStacks):
        argdata.append((function, parameter, subStacks[i]));    
    
    #run sequentially
    results = [];
    for i in range(nSubStacks):
        results.append(processSubStack(argdata[i]));
    
    #join the results
    results = join(results, subStacks = subStacks, **parameter);
    
    #write / or return 
    return io.writePoints(sink, results);












### Pickle does not like classes:

## sub stack information
#class SubStack(object):
#    """Class containing all info of a sub stack usefull for the image processing and result joining functions"""
#    
#    # sub stack id
#    stackId = None;
#    
#    # number of stacks
#    nStacks = None;
#    
#    # tuple of x,y,z range of this sub stack
#    z = all; 
#    x = all;
#    y = all;
#    
#    #original source
#    source = None;    
#    
#    # tuple of center point of the overlaping regions
#    zCenters = None;
#    
#    # tuple of z indices that would generate full image without overlaps
#    zCenterIndices = None;
#    
#    # tuple of z indices in the sub image as returned by readData that would generate full image without overlaps
#    zSubCenterIndices = None;
#    
#    
#    def __init__(slf, stackId = 0, nStacks = 1, source = None, x = all, y = all, z = all, zCenters = all, zCenterIndices = all):
#        slf.stackId = stackId;
#        slf.nStacks = nStacks;
#        slf.source = source; 
#        slf.x = x;
#        slf.y = y;
#        slf.z = z;
#        slf.zCenters = zCenters;
#        slf.zCenterIndices = zCenterIndices;
#        if not zCenterIndices is all and not z is all:
#            slf.zSubCenterIndices = (c - z[0] for c in zCenterIndices);
#        else:
#            slf.zSubCenterIndices = all;
#       
#    
#def printSubStackInfo(slf, out = sys.stdout):
#    out.write("Sub Stack: %d / %d\n" % (slf.stackId, slf.nStacks));
#    out.write("source:         %s\n" %       slf.source);
#    out.write("x,y,z:          %s, %s, %s\n" % (str(slf.x), str(slf.y), str(slf.z)));
#    out.write("zCenters:       %s\n" %       str(slf.zCenters));   
#    out.write("zCenterIndices: %s\n" %       str(slf.zCenterIndices));
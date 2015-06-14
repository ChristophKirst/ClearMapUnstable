# -*- coding: utf-8 -*-

import sys
self = sys.modules[__name__];

import math
import numpy

from multiprocessing import Pool

import iDISCO.IO.Imaris as io
from iDISCO.ImageProcessing.SpotDetection import detectCells;
from iDISCO.Utils.ProcessWriter import ProcessWriter;
from iDISCO.Utils.Timer import Timer;

#define the subroutine for the processing

def processSubStack(dsr):
    
    fn = dsr[0];
    segmentation = dsr[1];
    zrange  = dsr[2];
    iid = dsr[3];
    n = dsr[4];
    
    pw = ProcessWriter(iid);
    pw.write("processing chunk " + str(iid) + "/" + str(n));
    pw.write("file   = " + fn);
    pw.write("seg    = " + str(segmentation));
    pw.write("zrange = " + str(zrange)); 
    
    timer = Timer();
    
    img = readData(fn, zrange = zrange);
    
    pw.write(timer.elapsedTime(head = 'Reading data of size ' + str(img.shape)));
    timer.reset();    
    
    
    seg = segmentation(img, out = pw);
    
    pw.write(timer.elapsedTime(head = 'Segmenting data of size ' + str(img.shape)));  
    
    return seg;
    #return 0;
    #return numpy.array([[0,0,0.], [1,1,1]]);

    
def parallelProcessStack(fn, chunksizemax = 100, chunksizemin = 30, chunkoverlap = 15, processes = 2, segmentation = detectCells, zrange = all):
    """Parallel segmetation on a stack"""
    
    #determine z ranges
    if zrange == all: 
        f = io.openFile(fn);
        dataset = io.readData(f, resolution=0);
        nz = dataset.shape[2];    
        f.close();
        zrange = (0,nz);
    else:
        nz = zrange[1] - zrange[0];
    
    #calcualte chunk sizes
    chunksize = chunksizemax;
    nchunks = int(math.ceil((nz - chunksize) / (1. * (chunksize - chunkoverlap)) + 1));   
    chunksize = (nz + (nchunks-1) * chunkoverlap) / nchunks;
       
    #increase overlap if chunks to small
    chunksizemin = min(chunksizemin, chunkoverlap);
    if chunksize < chunksizemin:
        print "Warning: parallelProcessStack: optimal chunk size " + str(chunksize) + " smaller than minimum chunk size " + str(chunksizemin) + "!"; 
        chunksize = chunksizemin;
        chunkoverlap = math.ceil(chunksize - (nz - chunksize) / (nchunks -1));
        print "Warning: parallelProcessStack: setting chunk overlap to " + str(chunkoverlap) + "!";
           
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
            zhi = nz;
        
        zranges.append((int(zlo), int(zhi)));
        zcenters.append((zhiold - zlo) / 2. + zlo); 
        
    zcenters.append(nz);
    
    #adjust for the zrange
    zcenters = [x + zrange[0] for x in zcenters];
    zranges = [(x[0] + zrange[0], x[1] + zrange[0]) for x in zranges];
    
    #print zcenters
    #print zranges    
    
    
    argdata = []
    for i in range(nchunks):
        argdata.append((fn, segmentation, zranges[i], i, nchunks));
        
    print argdata
    
    # process in parallel
    pool = Pool(processes = processes);
    print processes
    
    resultsseg = pool.map(processSubStack, argdata);
    
    #print '=========== results';
    #print resultsseg;
    #print nchunks
    
    #print zranges
    #print zcenters
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
            
    
    return (numpy.concatenate(results), numpy.concatenate(resultsi))
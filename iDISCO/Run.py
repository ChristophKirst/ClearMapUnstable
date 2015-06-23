# -*- coding: utf-8 -*-
"""
Routines to run iDISCO alignment and cell detection

Created on Fri Jun 19 12:20:18 2015

@author: ckirst
"""

#import os

from iDISCO.Parameter import *
#from iDISCO.IO.IO import writePoints

import iDISCO.ImageProcessing.SpotDetection

try:
    import iDISCO.ImageProcessing.IlastikClassification
    haveIlastik = True;
except:
    haveIlastik = False;

from iDISCO.ImageProcessing.StackProcessing import parallelProcessStack, sequentiallyProcessStack

from iDISCO.Alignment.Elastix import transformPoints, alignData

from iDISCO.Utils.Timer import Timer

    

def runCellDetection(parameter):
    timer = Timer();
    
    pp = parameter.StackProcessing;
    ps = parameter.DataSource;
    
    # run segmentation
    if parameter.ImageProcessing.Method == "SpotDetection":
        
        detectCells = iDISCO.ImageProcessing.SpotDetection.detectCells;
        
        centers, intensities = parallelProcessStack(ps.ImageFile, x = ps.XRange, y = ps.YRange, z = ps.ZRange, 
                                                processes = pp.Processes, chunksizemax = pp.ChunkSizeMax, chunksizemin = pp.ChunkSizeMin, chunkoverlap = pp.ChunkOverlap, 
                                                optimizechunks = pp.OptimizeChunks, optimizechunksizeincrease = pp.OptimizeChunkSizeIncrease,
                                                segmentation = detectCells, parameter = parameter.ImageProcessing);        
        
    else:
        if haveIlastik:
            #ilastik does parallel processing so do sequential processing here
            detectCells = iDISCO.ImageProcessing.IlastikClassification.detectCells;
            
            centers, intensities = sequentiallyProcessStack(ps.ImageFile, x = ps.XRange, y = ps.YRange, z = ps.ZRange, 
                                                        chunksizemax = pp.ChunkSizeMax, chunksizemin = pp.ChunkSizeMin, chunkoverlap = pp.ChunkOverlap,
                                                        segmentation = detectCells, parameter = parameter.ImageProcessing);
            
        else:
            raise RuntimeError("No Ilastik installed use SpotDectection instead!");
        
    timer.printElapsedTime("Main");
    
    return centers, intensities;



def runCellCoordinateTransformation(parameter):
    
    cf = parameter.ImageProcessing.CellCoordinateFile;
    ap = parameter.Alignment;
    
    # downscale points to referenece image size
    print 'todo'    
    
    # transform points
    points = transformPoints(cf, transformparameterfile = None, parameter = ap, read = True, tmpfile = None, outdir = None, indices = False);
       
    # upscale ppints back to original size
    print 'todo'
    
    return points;
    
    
    
def runAlignment(parameter):
    
    fi = parameter.Alignmnet.FixedImage;
    mi = parameter.Alignment.MovingImage;
    
    alignData(mi, fi, parameter = parameter);
        
    return parameter.Alignment.AlignmentDirectory;
    
    
    



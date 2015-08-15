# -*- coding: utf-8 -*-
"""
Cell Detection Module

Notes:
    - main script to start the individual cell detection routines:
        - spot detection
        - illastik
        - ...
        
Created on Fri Aug 14 03:34:23 2015

@author: ckirst
"""

from iDISCO.IO import IO as io
import iDISCO.ImageProcessing.SpotDetection

try:
    import iDISCO.ImageProcessing.IlastikClassification
    haveIlastik = True;
except:
    haveIlastik = False;

from iDISCO.ImageProcessing.StackProcessing import parallelProcessStack, sequentiallyProcessStack

from iDISCO.Utils.Timer import Timer
    

def detectCells(source, method = **parameter):
    """Detect cells in data"""
    
    timer = Timer();
        
    # run segmentation
    if method == "SpotDetection":
        
        detectCells = iDISCO.ImageProcessing.SpotDetection.detectCells;
        
        centers, intensities = parallelProcessStack(source, x = ps.XRange, y = ps.YRange, z = ps.ZRange, 
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
    
    if not parameter.ImageProcessing.Parameter.ThresholdSave is None:
        iid = intensities >  parameter.ImageProcessing.Parameter.ThresholdSave;
        centers = centers[iid,:];

    if not parameter.ImageProcessing.PointFile is None:
        io.writePoints(parameter.ImageProcessing.PointFile, centers);
        
    if not parameter.ImageProcessing.IntensityFile is None:
        io.writePoints(parameter.ImageProcessing.IntensityFile, intensities);

    return centers, intensities;


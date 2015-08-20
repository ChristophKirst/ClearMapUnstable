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

import iDISCO.ImageProcessing.SpotDetection
import iDISCO.ImageProcessing.IlastikClassification

from iDISCO.ImageProcessing.StackProcessing import parallelProcessStack, sequentiallyProcessStack

from iDISCO.Utils.Timer import Timer
    

def detectCells(source, sink = None, method ="SpotDetection", processMethod = all, **parameter):
    """Detect cells in data"""
    timer = Timer();
        
    # run segmentation
    if method == "SpotDetection":
        detectCells = iDISCO.ImageProcessing.SpotDetection.detectCells;
    elif method == 'Ilastik':
        if iDISCO.ImageProcessing.IlastikClassification.Initialized:
            detectCells = iDISCO.ImageProcessing.IlastikClassification.detectCells;
            processMethod = 'sequential';  #ilastik does parallel processing so force sequential processing here
        else:
            raise RuntimeError("detectCells: Ilastik not initialized, fix in Settings.py or use SpotDectection method instead!");
    else:
        raise RuntimeError("detectCells: invalid method %s" % str(method));
    
    if processMethod == 'sequential':
        result = sequentiallyProcessStack(source, sink = sink, function = detectCells, **parameter);  
    elif processMethod is all or processMethod == 'parallel':
        result = parallelProcessStack(source, sink = sink, function = detectCells, **parameter);  
    else:
        raise RuntimeError("detectCells: invalid processMethod %s" % str(processMethod));
    
    timer.printElapsedTime("Total Cell Detection");
    
    return result;



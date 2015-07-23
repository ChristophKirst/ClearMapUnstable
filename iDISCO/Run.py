# -*- coding: utf-8 -*-
"""
Routines to run iDISCO alignment and cell detection

Created on Fri Jun 19 12:20:18 2015

@author: ckirst
"""

#import os

from iDISCO.Parameter import *
from iDISCO.IO import IO as io
import iDISCO.ImageProcessing.SpotDetection

try:
    import iDISCO.ImageProcessing.IlastikClassification
    haveIlastik = True;
except:
    haveIlastik = False;

from iDISCO.ImageProcessing.StackProcessing import parallelProcessStack, sequentiallyProcessStack

from iDISCO.Alignment.Elastix import transformPoints, alignData, initializeElastix

from iDISCO.Utils.Timer import Timer

from iDISCO.Alignment.Resampling import resampleData
    

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
    
    if not parameter.ImageProcessing.Parameter.ThresholdSave is None:
        iid = intensities >  parameter.ImageProcessing.Parameter.ThresholdSave;
        centers = centers[iid,:];

    if not parameter.ImageProcessing.CellCoordinateFile is None:
        io.writePoints(parameter.ImageProcessing.CellCoordinateFile, centers);
        
    if not parameter.ImageProcessing.CellIntensityFile is None:
        io.writePoints(parameter.ImageProcessing.CellIntensityFile, intensities);

    return centers, intensities;



def runInitializeElastix(parameter):

    ed = parameter.Alignment.ElastixDirectory;
    
    initializeElastix(ed);
    
    return ed;






def runCellCoordinateTransformation(parameter):
    
    cf = parameter.ImageProcessing.CellCoordinateFile;
    pa = parameter.Alignment;
    pr = parameter.Resampling;
    
    # downscale points to referenece image size
     
    
    
    # transform points
    points = transformPoints(cf, alignmentDirectory = ap.AlignmentDirectory, transformparameterfile = None, read = True, tmpfile = None, outdir = None, indices = False);
       
    # upscale ppints back to original size
    #print 'todo'
    
    return points;
    






    
def runAlignment(parameter):
    
    pa = parameter.Alignment;
    
    fi = pa.FixedImage;
    mi = pa.MovingImage;
    
    af = pa.AffineParameterFile;
    bf = pa.BSplineParameterFile;
    
    od = pa.AlignmentDirectory;
    
    alignData(mi, fi, af, bf, od);
        
    return od;
    
    
    
def runResampling(parameter):
    
    rp = parameter.Resampling;    
    
    resampleData(rp.DataFiles, outputFile = rp.ResampledFile, 
                 resolutionData = rp.ResolutionData, resolutionReference = rp.ResolutionReference,
                 processingDirectory = None, processes = rp.Processes, cleanup = True, orientation = rp.Orientation);
    
    return rp.ResampledFile;

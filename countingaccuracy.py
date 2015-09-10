# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 15:08:04 2015

@author: mtllab
"""


import os, numpy, math
import matplotlib.pyplot as pp
import iDISCO.Visualization.Plot as plt;
import iDISCO.Settings as settings
import iDISCO.IO.IO as io
import numpy as np
from iDISCO.Alignment.Resampling import resampleData;
from iDISCO.Alignment.Elastix import alignData, transformPoints
from iDISCO.ImageProcessing.CellDetection import detectCells
from iDISCO.Alignment.Resampling import resamplePoints, resamplePointsInverse
from iDISCO.Analysis.Label import countPointsInRegions
from iDISCO.Analysis.Voxelization import voxelize
from iDISCO.Analysis.Statistics import thresholdPoints
from iDISCO.Utils.ParameterTools import joinParameter
import iDISCO.Analysis.Label as lbl

######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/countingaccuracy';
   
cFosFile = '/home/mtllab/Documents/countingaccuracy/countingaccuracy.tif';
cFosFileRange = {'x' : all, 'y' : all, 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (815,1000), 'y' : (1078,1271), 'z' : (667,742)};


#Resolution of the Data (in um / pixel)
OriginalResolution = (4.0625, 4.0625, 3);


######################### Cell Detection Parameters using DoG

SpotDetectionParameter = {
    # background correctoin: None or (x,y) which is size of disk for gray scale opening
    "backgroundSize" : (15,15),
    
    # spot Detection via Difference of Gaussians (DoG) filter: (x,y,z) size
    "dogSize" : (7,7,11),
    
    # h of h-max transform
    "hMax" : 20,
    
    # intensity detection   
    "intensityMethod"  : 'Max',  #None -> intensity of pixel of center, alternatively string of numpy array method that returns a single number
    "intensitySize"    : (3,3,3),  # size of box in (x,y,z) to include in intensity determination
    
    # threshold for min intensity at center to be counted as cell, for saving ('None' will save everything )
    "threshold" : None,
      
    # write cell mask to disk (to check cell detection accuracy), if not None
    #"cellMaskFile" : os.path.join(BaseDirectory, 'cell_mask/cell_mask_Z\d{4}.ome.tif')
    "cellMaskFile" : None
    
    #some debug / quality check output
    #"verbose" : True,
    #"processMethod" : "sequential"  #  plotting during image processing only in sequential mode !
    };



############################ Config parameters

#Processes to use for Resampling
ResamplingParameter = { "processes": 12 };


#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 4,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "chunkOverlap" : 30,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : False,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
    "chunkOptimizationSize" : all
    };


# result files for cell coordinates (csv, vtk or ims)
SpotDetectionParameter["source"] = cFosFile;
SpotDetectionParameter = joinParameter(SpotDetectionParameter, cFosFileRange)

SpotDetectionParameter["sink"] = (os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy'));

ImageProcessingParameter = joinParameter(StackProcessingParameter, SpotDetectionParameter)


## Parameter to calculate density voxelization
VoxelizationParameter = {
    #Method to voxelize
    "voxelizationMethod" : 'Pixel', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized in pixels
    "voxelizationSize" : (1,1,1),  

    # Voxelization weigths (e/g intensities)
    "voxelizationWeights" : None
    };



######################################################
######################################################
######################################################

detectCells(**ImageProcessingParameter);
points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));


######################################################
######################################################
######################################################

label = io.readData('/home/mtllab/Documents/countingaccuracy/eliza.tif');
label = label.astype('int32');
humanpoint = label == 1;
check = [0,0,0];

for t in range(np.max(intensities[:,0]).astype(int)):
    pointT, intensitiesT = thresholdPoints(points, intensities, threshold = (t, 65000), row = (1,0));
    counts = lbl.countPointsInRegions(pointT, label);
    match = float(counts[1][1])/len(pointT);    
    check =  np.vstack([check,[t,float(match),len(pointT)]]);


######################################################
######################################################
######################################################



io.writePoints('/home/mtllab/Documents/countingaccuracy/nicocheckfiltered.csv', check);
pp.plot(check[:,0],check[:,1]);


######################################################
######################################################
######################################################

pointT, intensitiesT = thresholdPoints(points, intensities, threshold = (2200, 65000), row = (0,0));
io.writePoints('/home/mtllab/Documents/countingaccuracy/cells-2200_fast_unfiltered.npy', pointT)

######################################################
######################################################
######################################################


dataSource =  '/home/mtllab/Documents/countingaccuracy/countingaccuracy.tif'
pointSource= '/home/mtllab/Documents/countingaccuracy/cells-2200_fast_unfiltered.npy';
x = all;
y = all;
z = all;

vox = voxelize(pointSource, dataSource, **VoxelizationParameter);
io.writeData('/home/mtllab/Documents/countingaccuracy/countingaccuracy-2200_fast_unfiltered.tif', vox.astype('int32'));

######################################################
######################################################
######################################################
points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
label = io.readData('/home/mtllab/Documents/countingaccuracy/nico.nrrd');
label = label.astype('int32');
pointT, intensitiesT = thresholdPoints(points, intensities, threshold = (20, 65000), row = (1,0));
lblpoints = lbl.labelPoints(pointT, label);
pp.hist(intensitiesT[lblpoints==1][:,1], histtype = 'step', bins = 100, range=(10,200))
pp.hist(intensities[:,1], histtype = 'step', bins = 100, range=(10,200))


######################################################
######################################################
######################################################

labeleliza = io.readData('/home/mtllab/Documents/countingaccuracy/eliza.tif');
labeleliza = labeleliza.astype('int32');
labelnico = io.readData('/home/mtllab/Documents/countingaccuracy/nico.nrrd');
labelnico = labelnico.astype('int32');
elizapoint = labeleliza == 1;
nicopoint = labelnico == 1;

check = [0,0,0];

pointT, intensitiesT = thresholdPoints(points, intensities, threshold = (20, 65000), row = (1,0));
lblE = lbl.labelPoints(pointT, labeleliza);
lblN = lbl.labelPoints(pointT, labelnico);

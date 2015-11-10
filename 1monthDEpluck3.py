# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 19:50:50 2015

@author: mtllab
"""


import os, numpy, math

import iDISCO.Settings as settings
import iDISCO.IO.IO as io

from iDISCO.Alignment.Resampling import resampleData;
from iDISCO.Alignment.Elastix import alignData, transformPoints
from iDISCO.ImageProcessing.CellDetection import detectCells
from iDISCO.Alignment.Resampling import resamplePoints, resamplePointsInverse
from iDISCO.Analysis.Label import countPointsInRegions
from iDISCO.Analysis.Voxelization import voxelize
from iDISCO.Analysis.Statistics import thresholdPoints
from iDISCO.Utils.ParameterTools import joinParameter
from iDISCO.Analysis.Label import labelToName


######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/onemonthpluckDE/Pluck3';

cFosFile = '/home/mtllab/Documents/onemonthpluckDE/Pluck3/150914_0-8xs3-cfos20HFcont_19-27-28/19-27-28_0-8xs3-cfos20HFcont_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/onemonthpluckDE/Pluck3/150915_0-8xs3-autofluor_20-48-06/20-48-06_0-8xs3-autofluor_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (1000,1500), 'y' : (700,1200), 'z' : (500,1080)};


#Resolution of the Data (in um / pixel)
OriginalResolution = (4.0625, 4.0625, 3);

#Resolution of the Atlas (in um/ pixel)
AtlasResolution = (25, 25, 25);

#Orientation
FinalOrientation = (-1,2,3);

VoxelizationFile = os.path.join(BaseDirectory, 'points_voxelized.tif');



######################### Cell Detection Parameters using DoG

SpotDetectionParameter = {
    # background correctoin: None or (x,y) which is size of disk for gray scale opening
    "backgroundSize" : (15,15),
    
    # spot Detection via Difference of Gaussians (DoG) filter: (x,y,z) size
    #"dogSize" : (7, 7, 11),
    'dogSize' : None,
    
    # h of h-max transform
    "hMax" : None,
    
    # intensity detection   
    "intensityMethod"  : None, #'Max',  #None -> intensity of pixel of center, alternatively string of numpy array method that returns a single number
    "intensitySize"    : (3,3,3),  # size of box in (x,y,z) to include in intensity determination
    
    # threshold for min intensityNone at center to be counted as cell, for saving ('None' will save everything )
    "threshold" : None,
      
    # write cell mask to disk (to check cell detection accuracy), if not None
    #"cellMaskFile"   : os.path.join(BaseDirectory, 'cell_mask_new/cell_mask_Z\d{4}.ome.tif'),
    #"backgroundFile" : os.path.join(BaseDirectory, 'background/background_Z\d{4}.ome.tif'),
    #"hMaxFile"       : os.path.join(BaseDirectory, 'hmax/hmax_Z\d{4}.ome.tif'),
    #"dogFile" : os.path.join(BaseDirectory, 'dog/dog_Z\d{4}.ome.tif')
    #"cellMaskFile" : None
    
    #some debug / quality check output
    #"verbose" : True,
    #"processMethod" : "sequential"  #  plotting during image processing only in sequential mode !
    };



#################### Heat map generation

##Voxelization

## Parameter to calculate density voxelization
VoxelizationParameter = {
    #Method to voxelize
    "voxelizationMethod" : 'Spherical', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized in pixels
    "voxelizationSize" : (15,15,15),  

    # Voxelization weigths (e/g intensities)
    "voxelizationWeights" : None
    };


############################ Config parameters

#Processes to use for Resampling
ResamplingParameter = { "processes": 12 };

#Path to registration parameters and atlases
PathReg        = '/home/mtllab/Documents/warping';
AtlasFile      = os.path.join(PathReg, 'half_template_25_right.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right.tif');

#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 12,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 50,
    "chunkOverlap" : 32,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimization (True, False or all = automatic)
    "chunkOptimizationSize" : all
    };

ResolutionAffineCFosAutoFluo =  (16, 16, 16);




######################## Run Parameters

### Resample Fluorescent and CFos images
# Autofluorescent cFos resampling for aquisition correction
CorrectionResamplingParameterCfos = ResamplingParameter.copy();

CorrectionResamplingParameterCfos["source"] = cFosFile;
CorrectionResamplingParameterCfos["sink"]   = os.path.join(BaseDirectory, 'cfos_resampled.tif');
    
CorrectionResamplingParameterCfos["resolutionSource"] = OriginalResolution;
CorrectionResamplingParameterCfos["resolutionSink"]   = ResolutionAffineCFosAutoFluo;

CorrectionResamplingParameterCfos["orientation"] = FinalOrientation;
   
   
   
   #Files for Auto-fluorescence for acquisition movements correction
CorrectionResamplingParameterAutoFluo = CorrectionResamplingParameterCfos.copy();
CorrectionResamplingParameterAutoFluo["source"] = AutofluoFile;
CorrectionResamplingParameterAutoFluo["sink"]   = os.path.join(BaseDirectory, 'autofluo_for_cfos_resampled.tif');
   
#Files for Auto-fluorescence (Atlas Registration)
RegistrationResamplingParameter = CorrectionResamplingParameterCfos.copy();
RegistrationResamplingParameter["sink"]            =  os.path.join(BaseDirectory, 'autofluo_resampled.tif');
RegistrationResamplingParameter["resolutionSink"]  = AtlasResolution;
   

### Align cFos and Autofluo

CorrectionAlignmentParameter = {            
    #moving and reference images
    "movingImage" : os.path.join(BaseDirectory, 'autofluo_for_cfos_resampled.tif'),
    "fixedImage"  :  os.path.join(BaseDirectory, 'cfos_resampled.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(PathReg, 'Par0000affine_acquisition.txt'),
    "bSplineParameterFile" : None,
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_cfos_to_auto')
    }; 
  

### Align Autofluo and Atlas

#directory of the alignment result
RegistrationAlignmentParameter = CorrectionAlignmentParameter.copy();

RegistrationAlignmentParameter["resultDirectory"] = os.path.join(BaseDirectory, 'elastix_auto_to_atlas');
    
#moving and reference images
RegistrationAlignmentParameter["movingImage"]  = AtlasFile;
RegistrationAlignmentParameter["fixedImage"]   = os.path.join(BaseDirectory, 'autofluo_resampled.tif');

#elastix parameter files for alignment
RegistrationAlignmentParameter["affineParameterFile"]  = os.path.join(PathReg, 'Par0000affine.txt');
RegistrationAlignmentParameter["bSplineParameterFile"] = os.path.join(PathReg, 'Par0000bspline.txt');



# result files for cell coordinates (csv, vtk or ims)
SpotDetectionParameter["source"] = cFosFile;
SpotDetectionParameter = joinParameter(SpotDetectionParameter, cFosFileRange)

SpotDetectionParameter["sink"] = (os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy'));

ImageProcessingParameter = joinParameter(StackProcessingParameter, SpotDetectionParameter)


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to referenece image size

CorrectionResamplingPointsParameter = CorrectionResamplingParameterCfos.copy();
CorrectionResamplingPointsParameter["pointSource"] = os.path.join(BaseDirectory, 'cells.npy');
CorrectionResamplingPointsParameter["dataSizeSource"] = cFosFile;
CorrectionResamplingPointsParameter["pointSink"]  = None;

CorrectionResamplingPointsInverseParameter = CorrectionResamplingPointsParameter.copy();
CorrectionResamplingPointsInverseParameter["dataSizeSource"] = cFosFile;
CorrectionResamplingPointsInverseParameter["pointSink"]  = None;

## Transform points from corrected to registered
# downscale points to referenece image size
RegistrationResamplingPointParameter = RegistrationResamplingParameter.copy();
RegistrationResamplingPointParameter["dataSizeSource"] = cFosFile;
RegistrationResamplingPointParameter["pointSink"]  = None;


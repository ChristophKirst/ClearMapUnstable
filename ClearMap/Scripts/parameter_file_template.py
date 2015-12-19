# -*- coding: utf-8 -*-
"""
Example script to set up the parameters for the image processing pipeline
"""


import os, numpy, math

import ClearMap.Settings as settings
import ClearMap.IO as io

from ClearMap.Alignment.Resampling import resampleData;
from ClearMap.Alignment.Elastix import alignData, transformPoints
from ClearMap.ImageProcessing.CellDetection import detectCells
from ClearMap.Alignment.Resampling import resamplePoints, resamplePointsInverse
from ClearMap.Analysis.Label import countPointsInRegions
from ClearMap.Analysis.Voxelization import voxelize
from ClearMap.Analysis.Statistics import thresholdPoints
from ClearMap.Utils.ParameterTools import joinParameter
from ClearMap.Analysis.Label import labelToName


######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/whiskers/exploration/9';

cFosFile = '/home/mtllab/Documents/whiskers/exploration/9/150720_0_8xs3-cfos20HFcont_19-05-56/19-05-56_0_8xs3-cfos20HFcont_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/whiskers/exploration/9/150721_0_8xs3-autofluor_20-37-01/20-37-01_0_8xs3-autofluor_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
cFosFileRange = {'x' : (500,1000), 'y' : (900, 1400), 'z' : (500,1000)};

#Resolution of the Data (in um / pixel)
OriginalResolution = (4.0625, 4.0625, 3);

#Orientation
FinalOrientation = (1,2,3);

VoxelizationFile = os.path.join(BaseDirectory, 'points_voxelized.tif');


#Resolution of the Atlas (in um/ pixel)
AtlasResolution = (25, 25, 25);

#
ResolutionAffineCFosAutoFluo =  (16, 16, 16);

#Path to registration parameters and atlases
PathReg        = '/home/mtllab/Documents/warping';
AtlasFile      = os.path.join(PathReg, 'half_template_25_right.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right.tif');

#ImageProcessingMethod = "SpotDetection";

######################### Cell Detection Parameters using custom filters

correctIlluminationParameter = {
    "flatfield" : None,  # (str, True or None)  flat field intensities, if None d onot correct image for illumination, if True the 
    "background" : None, # (str, None or array) background image as file name or array, if None background is assumed to be zero
    "scaling" :  "Mean", # (str or None)        scale the corrected result by this factor, if 'max'/'mean' scale to keep max/mean invariant
    "save" : None,       # (str or None)        save the corrected image to file
    "verbose" : True    # (bool or int)        print / plot information about this step 
}

removeBackgroundParameter = {
    "size" : (7,7),  # size for the structure element of the morphological opening
    "save" : None,     # file name to save result of this operation
    "verbose" : True  # print / plot information about this step       
}


filterDoGParameter = {
    "size" :  None,        # (tuple or None)      size for the DoG filter if None, do not correct for any background
    "sigma" : None,        # (tuple or None)      std of outer Guassian, if None autmatically determined from size
    "sigma2": None,        # (tuple or None)      std of inner Guassian, if None autmatically determined from size
    "save"  : None,        # (str or None)        file name to save result of this operation if None dont save to file 
    "verbose" : True      # (bool or int)        print / plot information about this step
}

findExtendedMaximaParameter = {
    "hMax" : None,            # (float or None)     h parameter for the initial h-Max transform, if None, do not perform a h-max transform
    "size" : 5,             # (tuple)             size for the structure element for the local maxima filter
    "threshold" : 0,        # (float or None)     include only maxima larger than a threshold, if None keep all localmaxima
    "save"  : None,         # (str or None)       file name to save result of this operation if None dont save to file 
    "verbose" : True       # (bool or int)       print / plot information about this step
}

findIntensityParameter = {
    "method" : 'Max',       # (str, func, None)   method to use to determine intensity (e.g. "Max" or "Mean") if None take intensities at the given pixels
    "size" :  (3,3,3)       # (tuple)             size of the box on which to perform the *method*
}

detectCellShapeParameter = {
    "threshold" : 700,     # (float or None)      threshold to determine mask, pixel below this are background if None no mask is generated
    "save"  : None,        # (str or None)        file name to save result of this operation if None dont save to file 
    "verbose" : True      # (bool or int)        print / plot information about this step if None take intensities at the given pixels
}

## Paramters for cell detection using spot detection algorithm 
detectSpotsParameter = {
    "correctIlluminationParameter" : correctIlluminationParameter,
    "removeBackgroundParameter"    : removeBackgroundParameter,
    "filterDoGParameter"           : filterDoGParameter,
    "findExtendedMaximaParameter"  : findExtendedMaximaParameter,
    "findIntensityParameter"       : findIntensityParameter,
    "detectCellShapeParameter"     : detectCellShapeParameter
}





#################### Heat map generation

##Voxelization

## Parameter to calculate density voxelization
voxelizeParameter = {
    #Method to voxelize
    "method" : 'Spherical', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized in pixels
    "size" : (15,15,15),  

    # Voxelization weigths (e/g intensities)
    "weights" : None
};





############################ Config parameters

#Processes to use for Resampling
ResamplingParameter = {
    "processes": 12 
};


#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 6,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 50,
    "chunkOverlap" : 32,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimization (True, False or all = automatic)
    "chunkOptimizationSize" : all,
   
    "processMethod" : "parallel"   
   };






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
    "fixedImage"  : os.path.join(BaseDirectory, 'cfos_resampled.tif'),
    
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
SpotDetectionParameter = {
    "source" : cFosFile,
    "sink"   : (os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')),
    "detectSpotsParameter" : detectSpotsParameter
};
SpotDetectionParameter = joinParameter(SpotDetectionParameter, cFosFileRange)

ImageProcessingParameter = joinParameter(StackProcessingParameter, SpotDetectionParameter);

FilteredCellsFile = (os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy'));

TransformedCellsFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')

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


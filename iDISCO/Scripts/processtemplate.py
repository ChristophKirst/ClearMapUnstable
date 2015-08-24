
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 14:58:07 2015

@author: mtllab
"""


import os, numpy, math

import iDISCO.Settings as settings
import iDISCO.IO.IO as io

from iDISCO.Alignment.Resampling import resampleData;
from iDISCO.Alignment.Elastix import alignData
from iDISCO.ImageProcessing.CellDetection import detectCells
from iDISCO.Alignment.Resampling import resamplePoints, resamplePointsInverse
from iDISCO.Analysis.Label import countPointsInRegions

from iDISCO.Utils.ParameterTools import joinParameter

######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/whiskers/nonshaved/150620-6R';
   
cFosFile      = '/home/mtllab/Documents/Haloperidol/1267/150818_0_8X-cfos_17-14-34/17-14-34_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};

AutofluoFile  = '/home/mtllab/Documents/Haloperidol/1267/150818_0_8X-autofluo_18-15-14/18-15-14_0_8X-autofluo_UltraII_C00_xyz-Table Z\d{4}.ome.tif';

#Resolution of the Data (in um / pixel)
OriginalResolution = (4.0625, 4.0625, 3);

#Resolution of the Atlas (in um/ pixel)
AtlasResolution = (25, 25, 25);

#Orientation
FinalOrientation = (1,2,3);

VoxelizationFile = os.path.join(BaseDirectory, 'points_voxelized.tif');


######################### Cell Detection Parameters using DoG

SpotDetectionParameter = {
    # background correctoin: None or (x,y) which is size of disk for gray scale opening
    "backgroundSize" : (15,15),
    
    # spot Detection via Difference of Gaussians (DoG) filter: (x,y,z) size
    "dogSize" : (7, 7, 11),
    
    # h of h-max transform
    "hMax" : 20,
    
    # intensity detection   
    "intensityMethod"  : None,  #None -> intensity of pixel of center, alternatively string of numpy array method that returns a single number
    #"intensitySize"    : (3,3,3),  # size of box in (x,y,z) to include in intensity determination
    
    # threshold for min intensity at center to be counted as cell (should be similar to the h max)
    "threshold" : 20,
        
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
       
    # Define bounds of the volume to be voxelized
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
    "processes" : 5,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "cChunkOverlap" : 15,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
    "chunkOptimizationSize" : all
    };

ResolutionAffineCFosAutoFluo =  (16, 16, 16);



######################## RUN


### Resample Fluorescent and CFos images

# Autofluorescent cFos resampling for aquisition correction
CorrectionResamplingParameter = ResamplingParameter;

CorrectionResamplingParameter["source"] = cFosFile;
CorrectionResamplingParameter["sink"]   = os.path.join(BaseDirectory, 'cfos_resampled.tif');
    
CorrectionResamplingParameter["resolutionSource"] = OriginalResolution;
CorrectionResamplingParameter["resolutionSink"]   = ResolutionAffineCFosAutoFluo;

CorrectionResamplingParameter["orientation"] = FinalOrientation;
   
resampleData(**CorrectionResamplingParameter);


#Files for Auto-fluorescence for acquisition movements correction
CorrectionResamplingParameter["source"] = AutofluoFile;
CorrectionResamplingParameter["sink"]   = os.path.join(BaseDirectory, 'autofluo_for_cfos_resampled.tif');

resampleData(**CorrectionResamplingParameter);


#Files for Auto-fluorescence (Atlas Registration)
RegistrationResamplingParameter = CorrectionResamplingParameter;
RegistrationResamplingParameter["sink"]            =  os.path.join(BaseDirectory, 'autofluo_resampled.tif');
RegistrationResamplingParameter["resolutionSink"]  = AtlasResolution;

resampleData(**RegistrationResamplingParameter);



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

resultDirectory  = alignData(**CorrectionAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory


# Align Autofluo and Atlas

#directory of the alignment result
RegistrationAlignmentParameter = CorrectionAlignmentParameter;

RegistrationAlignmentParameter["resultDirectory"] = os.path.join(BaseDirectory, 'elastix_auto_to_atlas');
    
#moving and reference images
RegistrationAlignmentParameter["movingImage"]  = AtlasFile;
RegistrationAlignmentParameter["fixedImage"]   = os.path.join(BaseDirectory, 'autofluo_resampled.tif');

#elastix parameter files for alignment
RegistrationAlignmentParameter["affineParameterFile"]  = os.path.join(PathReg, 'Par0000affine.txt');
RegistrationAlignmentParameter["bSplineParameterFile"] = os.path.join(PathReg, 'Par0000bspline.txt');

resultDirectory  = alignData(**RegistrationAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory



### Detect Cells

## Cell Detection

# result files for cell coordinates (csv, vtk or ims)
SpotDetectionParameter["source"] = cFosFile;
SpotDetectionParameter = joinParameter(SpotDetectionParameter, cFosFileRange)

SpotDetectionParameter["sink"] = (os.path.join(BaseDirectory, 'cells.csv'),  os.path.join(BaseDirectory,  'intensities.csv'));

allParameter = joinParameter(StackProcessingParameter, SpotDetectionParameter)

detectCells(**allParameter);


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to referenece image size
points = resamplePoints(SpotDetectionParameter["sink"][0], cFosFile, **CorrectionResamplingParameter);
    
# transform points
points = resamplePoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

# resmaple back to original
points = resamplePointsInverse(points, cFosFile, **CorrectionResamplingParameter);


## Transform points from corrected to registered
# downscale points to referenece image size
points = resamplePoints(points, cFosFile, **RegistrationResamplingParameter);
    
# transform points
points = resamplePoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);


TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.csv')

io.writePoints(TransformedCellFile, points);



#################### Heat map generation


#points = io.readPoints(TransformedCellFile);
vox = voxelize(points, AtlasFile, VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));


##################### Label points by anatomical regions

#points = io.readPoints(TransformedCellFile);
counts = countPointsInRegions(points, AnnotationFile);
io.writePoints(os.path.join(BaseDirectory,'Annotated_Points.csv'), counts);

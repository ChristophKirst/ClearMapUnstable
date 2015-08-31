# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 20:23:41 2015

@author: mtllab
"""


# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 14:58:07 2015

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


######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/parental/1Galcre';
   
cFosFile = '/home/mtllab/Documents/parental/1Galcre/150618_0_8X-stainings_13-43-54/13-43-54_0_8X-stainings_UltraII_C00_UltraII Filter0000_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/parental/1Galcre/150618_0_8X-autofluo_19-36-40/19-36-40_0_8X-autofluo_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (815,1000), 'y' : (1078,1271), 'z' : (667,742)};


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

      
# Threshold for the points to be considered in the analysis
minthreshold = 10; #remove points detected in the background, based on the filtered intensities
maxthreshold = 40000; #remove staining artefacts (bright antibody clumps), based on the non-filtered original intensities 

     

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
AtlasFile      = os.path.join(PathReg, 'half_template_25_right_fullWD.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right_fullWD.tif');

#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 4,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "chunkOverlap" : 30,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
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
CorrectionResamplingPointsParameter["source"] = os.path.join(BaseDirectory, 'cells.npy');
CorrectionResamplingPointsParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsParameter["sink"]  = None;

CorrectionResamplingPointsInverseParameter = CorrectionResamplingPointsParameter.copy();
CorrectionResamplingPointsInverseParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsInverseParameter["sink"]  = None;

## Transform points from corrected to registered
# downscale points to referenece image size
RegistrationResamplingPointParameter = RegistrationResamplingParameter.copy();
RegistrationResamplingPointParameter["dataSize"] = cFosFile;
RegistrationResamplingPointParameter["sink"]  = None;








############################################################################### 
# RUN
###############################################################################   


### Resample Fluorescent and CFos images
   
resampleData(**CorrectionResamplingParameterCfos);

resampleData(**CorrectionResamplingParameterAutoFluo);

resampleData(**RegistrationResamplingParameter);


### Align cFos and Autofluo

resultDirectory  = alignData(**CorrectionAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory


### Align Autofluo and Atlas

resultDirectory  = alignData(**RegistrationAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory



### Detect Cells

## Cell Detection

detectCells(**ImageProcessingParameter);

## Save spot retained as legitimate cells, based on intensities:

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minthreshold, maxthreshold), row = (1,0));

io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to reference image size

points = io.readPoints(CorrectionResamplingPointsParameter["source"]);


points = resamplePoints(**CorrectionResamplingPointsParameter);
    
# transform points
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

# resmaple back to original
CorrectionResamplingPointsInverseParameter["source"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);

## Transform points from corrected to registered
# downscale points to reference image size
RegistrationResamplingPointParameter["source"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
    
# transform points
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')

io.writePoints(TransformedCellFile, points);



#################### Heat map generation

#points = io.readPoints(TransformedCellFile);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

##################### Label points by anatomical regions

#points = io.readPoints(TransformedCellFile);
counts = countPointsInRegions(points, AnnotationFile);
io.writePoints(os.path.join(BaseDirectory,'Annotated_Points.csv'), counts);

#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################



######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/parental/2Galcre';
   
cFosFile = '/home/mtllab/Documents/parental/2Galcre/150729_0_8X-cfos_17-53-01/17-53-01_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/parental/2Galcre/150729_0_8X-autofluo_19-54-08/19-54-08_0_8X-autofluo_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (815,1000), 'y' : (1078,1271), 'z' : (667,742)};


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

      
# Threshold for the points to be considered in the analysis
minthreshold = 10; #remove points detected in the background, based on the filtered intensities
maxthreshold = 40000; #remove staining artefacts (bright antibody clumps), based on the non-filtered original intensities 

     

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
AtlasFile      = os.path.join(PathReg, 'half_template_25_right_fullWD.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right_fullWD.tif');

#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 4,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "chunkOverlap" : 30,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
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
CorrectionResamplingPointsParameter["source"] = os.path.join(BaseDirectory, 'cells.npy');
CorrectionResamplingPointsParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsParameter["sink"]  = None;

CorrectionResamplingPointsInverseParameter = CorrectionResamplingPointsParameter.copy();
CorrectionResamplingPointsInverseParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsInverseParameter["sink"]  = None;

## Transform points from corrected to registered
# downscale points to referenece image size
RegistrationResamplingPointParameter = RegistrationResamplingParameter.copy();
RegistrationResamplingPointParameter["dataSize"] = cFosFile;
RegistrationResamplingPointParameter["sink"]  = None;








############################################################################### 
# RUN
###############################################################################   


### Resample Fluorescent and CFos images
   
resampleData(**CorrectionResamplingParameterCfos);

resampleData(**CorrectionResamplingParameterAutoFluo);

resampleData(**RegistrationResamplingParameter);


### Align cFos and Autofluo

resultDirectory  = alignData(**CorrectionAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory


### Align Autofluo and Atlas

resultDirectory  = alignData(**RegistrationAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory



### Detect Cells

## Cell Detection

detectCells(**ImageProcessingParameter);

## Save spot retained as legitimate cells, based on intensities:

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minthreshold, maxthreshold), row = (1,0));

io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to reference image size

points = io.readPoints(CorrectionResamplingPointsParameter["source"]);


points = resamplePoints(**CorrectionResamplingPointsParameter);
    
# transform points
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

# resmaple back to original
CorrectionResamplingPointsInverseParameter["source"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);

## Transform points from corrected to registered
# downscale points to reference image size
RegistrationResamplingPointParameter["source"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
    
# transform points
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')

io.writePoints(TransformedCellFile, points);



#################### Heat map generation

#points = io.readPoints(TransformedCellFile);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

##################### Label points by anatomical regions

#points = io.readPoints(TransformedCellFile);
counts = countPointsInRegions(points, AnnotationFile);
io.writePoints(os.path.join(BaseDirectory,'Annotated_Points.csv'), counts);

#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################


######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/parental/3Galcre';
   
cFosFile = '/home/mtllab/Documents/parental/3Galcre/150728_0_8X-cfos_18-24-44/18-24-44_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/parental/3Galcre/150729_0_8X-autofluo_20-28-04/20-28-04_0_8X-autofluo_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (815,1000), 'y' : (1078,1271), 'z' : (667,742)};


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

      
# Threshold for the points to be considered in the analysis
minthreshold = 10; #remove points detected in the background, based on the filtered intensities
maxthreshold = 40000; #remove staining artefacts (bright antibody clumps), based on the non-filtered original intensities 

     

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
AtlasFile      = os.path.join(PathReg, 'half_template_25_right_fullWD.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right_fullWD.tif');

#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 4,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "chunkOverlap" : 30,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
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
CorrectionResamplingPointsParameter["source"] = os.path.join(BaseDirectory, 'cells.npy');
CorrectionResamplingPointsParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsParameter["sink"]  = None;

CorrectionResamplingPointsInverseParameter = CorrectionResamplingPointsParameter.copy();
CorrectionResamplingPointsInverseParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsInverseParameter["sink"]  = None;

## Transform points from corrected to registered
# downscale points to referenece image size
RegistrationResamplingPointParameter = RegistrationResamplingParameter.copy();
RegistrationResamplingPointParameter["dataSize"] = cFosFile;
RegistrationResamplingPointParameter["sink"]  = None;








############################################################################### 
# RUN
###############################################################################   


### Resample Fluorescent and CFos images
   
resampleData(**CorrectionResamplingParameterCfos);

resampleData(**CorrectionResamplingParameterAutoFluo);

resampleData(**RegistrationResamplingParameter);


### Align cFos and Autofluo

resultDirectory  = alignData(**CorrectionAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory


### Align Autofluo and Atlas

resultDirectory  = alignData(**RegistrationAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory



### Detect Cells

## Cell Detection

detectCells(**ImageProcessingParameter);

## Save spot retained as legitimate cells, based on intensities:

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minthreshold, maxthreshold), row = (1,0));

io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to reference image size

points = io.readPoints(CorrectionResamplingPointsParameter["source"]);


points = resamplePoints(**CorrectionResamplingPointsParameter);
    
# transform points
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

# resmaple back to original
CorrectionResamplingPointsInverseParameter["source"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);

## Transform points from corrected to registered
# downscale points to reference image size
RegistrationResamplingPointParameter["source"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
    
# transform points
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')

io.writePoints(TransformedCellFile, points);



#################### Heat map generation

#points = io.readPoints(TransformedCellFile);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

##################### Label points by anatomical regions

#points = io.readPoints(TransformedCellFile);
counts = countPointsInRegions(points, AnnotationFile);
io.writePoints(os.path.join(BaseDirectory,'Annotated_Points.csv'), counts);

#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################


######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/parental/4Galcre';
   
cFosFile = '/home/mtllab/Documents/parental/4Galcre/150729_0_8X-cfos_10-17-45/10-17-45_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/parental/4Galcre/150729_0_8X-autofluo_12-28-20/12-28-20_0_8X-autofluo_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (815,1000), 'y' : (1078,1271), 'z' : (667,742)};


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

      
# Threshold for the points to be considered in the analysis
minthreshold = 10; #remove points detected in the background, based on the filtered intensities
maxthreshold = 40000; #remove staining artefacts (bright antibody clumps), based on the non-filtered original intensities 

     

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
AtlasFile      = os.path.join(PathReg, 'half_template_25_right_fullWD.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right_fullWD.tif');

#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 4,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "chunkOverlap" : 30,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
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
CorrectionResamplingPointsParameter["source"] = os.path.join(BaseDirectory, 'cells.npy');
CorrectionResamplingPointsParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsParameter["sink"]  = None;

CorrectionResamplingPointsInverseParameter = CorrectionResamplingPointsParameter.copy();
CorrectionResamplingPointsInverseParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsInverseParameter["sink"]  = None;

## Transform points from corrected to registered
# downscale points to referenece image size
RegistrationResamplingPointParameter = RegistrationResamplingParameter.copy();
RegistrationResamplingPointParameter["dataSize"] = cFosFile;
RegistrationResamplingPointParameter["sink"]  = None;








############################################################################### 
# RUN
###############################################################################   


### Resample Fluorescent and CFos images
   
resampleData(**CorrectionResamplingParameterCfos);

resampleData(**CorrectionResamplingParameterAutoFluo);

resampleData(**RegistrationResamplingParameter);


### Align cFos and Autofluo

resultDirectory  = alignData(**CorrectionAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory


### Align Autofluo and Atlas

resultDirectory  = alignData(**RegistrationAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory



### Detect Cells

## Cell Detection

detectCells(**ImageProcessingParameter);

## Save spot retained as legitimate cells, based on intensities:

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minthreshold, maxthreshold), row = (1,0));

io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to reference image size

points = io.readPoints(CorrectionResamplingPointsParameter["source"]);


points = resamplePoints(**CorrectionResamplingPointsParameter);
    
# transform points
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

# resmaple back to original
CorrectionResamplingPointsInverseParameter["source"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);

## Transform points from corrected to registered
# downscale points to reference image size
RegistrationResamplingPointParameter["source"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
    
# transform points
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')

io.writePoints(TransformedCellFile, points);



#################### Heat map generation

#points = io.readPoints(TransformedCellFile);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

##################### Label points by anatomical regions

#points = io.readPoints(TransformedCellFile);
counts = countPointsInRegions(points, AnnotationFile);
io.writePoints(os.path.join(BaseDirectory,'Annotated_Points.csv'), counts);

#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################


######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/parental/5Galcre';
   
cFosFile = '/home/mtllab/Documents/parental/5Galcre/150701_0_8X-stainings_17-26-26/17-26-26_0_8X-stainings_UltraII_C00_UltraII Filter0000_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/parental/5Galcre/150702_0_8X-autofluo_23-23-45/23-23-45_0_8X-autofluo_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (815,1000), 'y' : (1078,1271), 'z' : (667,742)};


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

      
# Threshold for the points to be considered in the analysis
minthreshold = 10; #remove points detected in the background, based on the filtered intensities
maxthreshold = 40000; #remove staining artefacts (bright antibody clumps), based on the non-filtered original intensities 

     

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
AtlasFile      = os.path.join(PathReg, 'half_template_25_right_fullWD.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right_fullWD.tif');

#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 4,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "chunkOverlap" : 30,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
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
CorrectionResamplingPointsParameter["source"] = os.path.join(BaseDirectory, 'cells.npy');
CorrectionResamplingPointsParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsParameter["sink"]  = None;

CorrectionResamplingPointsInverseParameter = CorrectionResamplingPointsParameter.copy();
CorrectionResamplingPointsInverseParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsInverseParameter["sink"]  = None;

## Transform points from corrected to registered
# downscale points to referenece image size
RegistrationResamplingPointParameter = RegistrationResamplingParameter.copy();
RegistrationResamplingPointParameter["dataSize"] = cFosFile;
RegistrationResamplingPointParameter["sink"]  = None;








############################################################################### 
# RUN
###############################################################################   


### Resample Fluorescent and CFos images
   
resampleData(**CorrectionResamplingParameterCfos);

resampleData(**CorrectionResamplingParameterAutoFluo);

resampleData(**RegistrationResamplingParameter);


### Align cFos and Autofluo

resultDirectory  = alignData(**CorrectionAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory


### Align Autofluo and Atlas

resultDirectory  = alignData(**RegistrationAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory



### Detect Cells

## Cell Detection

detectCells(**ImageProcessingParameter);

## Save spot retained as legitimate cells, based on intensities:

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minthreshold, maxthreshold), row = (1,0));

io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to reference image size

points = io.readPoints(CorrectionResamplingPointsParameter["source"]);


points = resamplePoints(**CorrectionResamplingPointsParameter);
    
# transform points
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

# resmaple back to original
CorrectionResamplingPointsInverseParameter["source"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);

## Transform points from corrected to registered
# downscale points to reference image size
RegistrationResamplingPointParameter["source"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
    
# transform points
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')

io.writePoints(TransformedCellFile, points);



#################### Heat map generation

#points = io.readPoints(TransformedCellFile);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

##################### Label points by anatomical regions

#points = io.readPoints(TransformedCellFile);
counts = countPointsInRegions(points, AnnotationFile);
io.writePoints(os.path.join(BaseDirectory,'Annotated_Points.csv'), counts);

#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################


######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/parental/id2Ucn3cre';
   
cFosFile = '/home/mtllab/Documents/parental/id2Ucn3cre/150727_0_8X-cfos_16-42-20/16-42-20_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/parental/id2Ucn3cre/150727_0_8X-autofluo_18-55-59/18-55-59_0_8X-autofluo_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (815,1000), 'y' : (1078,1271), 'z' : (667,742)};


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
    "dogSize" : (7, 7, 11),
    
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

      
# Threshold for the points to be considered in the analysis
minthreshold = 10; #remove points detected in the background, based on the filtered intensities
maxthreshold = 40000; #remove staining artefacts (bright antibody clumps), based on the non-filtered original intensities 

     

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
AtlasFile      = os.path.join(PathReg, 'half_template_25_right_fullWD.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right_fullWD.tif');

#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 4,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "chunkOverlap" : 30,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
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
CorrectionResamplingPointsParameter["source"] = os.path.join(BaseDirectory, 'cells.npy');
CorrectionResamplingPointsParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsParameter["sink"]  = None;

CorrectionResamplingPointsInverseParameter = CorrectionResamplingPointsParameter.copy();
CorrectionResamplingPointsInverseParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsInverseParameter["sink"]  = None;

## Transform points from corrected to registered
# downscale points to referenece image size
RegistrationResamplingPointParameter = RegistrationResamplingParameter.copy();
RegistrationResamplingPointParameter["dataSize"] = cFosFile;
RegistrationResamplingPointParameter["sink"]  = None;








############################################################################### 
# RUN
###############################################################################   


### Resample Fluorescent and CFos images
   
resampleData(**CorrectionResamplingParameterCfos);

resampleData(**CorrectionResamplingParameterAutoFluo);

resampleData(**RegistrationResamplingParameter);


### Align cFos and Autofluo

resultDirectory  = alignData(**CorrectionAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory


### Align Autofluo and Atlas

resultDirectory  = alignData(**RegistrationAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory



### Detect Cells

## Cell Detection

detectCells(**ImageProcessingParameter);

## Save spot retained as legitimate cells, based on intensities:

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minthreshold, maxthreshold), row = (1,0));

io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to reference image size

points = io.readPoints(CorrectionResamplingPointsParameter["source"]);


points = resamplePoints(**CorrectionResamplingPointsParameter);
    
# transform points
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

# resmaple back to original
CorrectionResamplingPointsInverseParameter["source"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);

## Transform points from corrected to registered
# downscale points to reference image size
RegistrationResamplingPointParameter["source"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
    
# transform points
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')

io.writePoints(TransformedCellFile, points);



#################### Heat map generation

#points = io.readPoints(TransformedCellFile);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

##################### Label points by anatomical regions

#points = io.readPoints(TransformedCellFile);
counts = countPointsInRegions(points, AnnotationFile);
io.writePoints(os.path.join(BaseDirectory,'Annotated_Points.csv'), counts);

#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################


######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/parental/id3Ucn3cre';
   
cFosFile = '/home/mtllab/Documents/parental/id3Ucn3cre/150728_0_8X-cfos_22-02-44/22-02-44_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/parental/id3Ucn3cre/150728_0_8X-autofluo_00-15-17/00-15-17_0_8X-autofluo_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (815,1000), 'y' : (1078,1271), 'z' : (667,742)};


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
    "dogSize" : (7, 7, 11),
    
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

      
# Threshold for the points to be considered in the analysis
minthreshold = 10; #remove points detected in the background, based on the filtered intensities
maxthreshold = 40000; #remove staining artefacts (bright antibody clumps), based on the non-filtered original intensities 

     

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
AtlasFile      = os.path.join(PathReg, 'half_template_25_right_fullWD.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right_fullWD.tif');

#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 4,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "chunkOverlap" : 30,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
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
CorrectionResamplingPointsParameter["source"] = os.path.join(BaseDirectory, 'cells.npy');
CorrectionResamplingPointsParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsParameter["sink"]  = None;

CorrectionResamplingPointsInverseParameter = CorrectionResamplingPointsParameter.copy();
CorrectionResamplingPointsInverseParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsInverseParameter["sink"]  = None;

## Transform points from corrected to registered
# downscale points to referenece image size
RegistrationResamplingPointParameter = RegistrationResamplingParameter.copy();
RegistrationResamplingPointParameter["dataSize"] = cFosFile;
RegistrationResamplingPointParameter["sink"]  = None;








############################################################################### 
# RUN
###############################################################################   


### Resample Fluorescent and CFos images
   
resampleData(**CorrectionResamplingParameterCfos);

resampleData(**CorrectionResamplingParameterAutoFluo);

resampleData(**RegistrationResamplingParameter);


### Align cFos and Autofluo

resultDirectory  = alignData(**CorrectionAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory


### Align Autofluo and Atlas

resultDirectory  = alignData(**RegistrationAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory



### Detect Cells

## Cell Detection

detectCells(**ImageProcessingParameter);

## Save spot retained as legitimate cells, based on intensities:

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minthreshold, maxthreshold), row = (1,0));

io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to reference image size

points = io.readPoints(CorrectionResamplingPointsParameter["source"]);


points = resamplePoints(**CorrectionResamplingPointsParameter);
    
# transform points
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

# resmaple back to original
CorrectionResamplingPointsInverseParameter["source"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);

## Transform points from corrected to registered
# downscale points to reference image size
RegistrationResamplingPointParameter["source"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
    
# transform points
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')

io.writePoints(TransformedCellFile, points);



#################### Heat map generation

#points = io.readPoints(TransformedCellFile);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

##################### Label points by anatomical regions

#points = io.readPoints(TransformedCellFile);
counts = countPointsInRegions(points, AnnotationFile);
io.writePoints(os.path.join(BaseDirectory,'Annotated_Points.csv'), counts);

#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################


######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/parental/id4Ucn3cre';
   
cFosFile = '/home/mtllab/Documents/parental/id4Ucn3cre/150727_0_8X-cfos_19-23-19/19-23-19_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/parental/id4Ucn3cre/150728_0_8X-autofluo_21-36-19/21-36-19_0_8X-autofluo_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (815,1000), 'y' : (1078,1271), 'z' : (667,742)};


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
    "dogSize" : (7, 7, 11),
    
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

      
# Threshold for the points to be considered in the analysis
minthreshold = 10; #remove points detected in the background, based on the filtered intensities
maxthreshold = 40000; #remove staining artefacts (bright antibody clumps), based on the non-filtered original intensities 

     

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
AtlasFile      = os.path.join(PathReg, 'half_template_25_right_fullWD.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right_fullWD.tif');

#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 4,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "chunkOverlap" : 30,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
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
CorrectionResamplingPointsParameter["source"] = os.path.join(BaseDirectory, 'cells.npy');
CorrectionResamplingPointsParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsParameter["sink"]  = None;

CorrectionResamplingPointsInverseParameter = CorrectionResamplingPointsParameter.copy();
CorrectionResamplingPointsInverseParameter["dataSize"] = cFosFile;
CorrectionResamplingPointsInverseParameter["sink"]  = None;

## Transform points from corrected to registered
# downscale points to referenece image size
RegistrationResamplingPointParameter = RegistrationResamplingParameter.copy();
RegistrationResamplingPointParameter["dataSize"] = cFosFile;
RegistrationResamplingPointParameter["sink"]  = None;








############################################################################### 
# RUN
###############################################################################   


### Resample Fluorescent and CFos images
   
resampleData(**CorrectionResamplingParameterCfos);

resampleData(**CorrectionResamplingParameterAutoFluo);

resampleData(**RegistrationResamplingParameter);


### Align cFos and Autofluo

resultDirectory  = alignData(**CorrectionAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory


### Align Autofluo and Atlas

resultDirectory  = alignData(**RegistrationAlignmentParameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory



### Detect Cells

## Cell Detection

detectCells(**ImageProcessingParameter);

## Save spot retained as legitimate cells, based on intensities:

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minthreshold, maxthreshold), row = (1,0));

io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to reference image size

points = io.readPoints(CorrectionResamplingPointsParameter["source"]);


points = resamplePoints(**CorrectionResamplingPointsParameter);
    
# transform points
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

# resmaple back to original
CorrectionResamplingPointsInverseParameter["source"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);

## Transform points from corrected to registered
# downscale points to reference image size
RegistrationResamplingPointParameter["source"] = points;
points = resamplePoints(**RegistrationResamplingPointParameter);
    
# transform points
points = transformPoints(points, transformDirectory = RegistrationAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

TransformedCellFile = os.path.join(BaseDirectory, 'cells_transformed_to_Atlas.npy')

io.writePoints(TransformedCellFile, points);



#################### Heat map generation

#points = io.readPoints(TransformedCellFile);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap.tif'), vox.astype('int32'));

VoxelizationParameter["voxelizationWeights"] = intensities[:,0].astype(float);
vox = voxelize(points, AtlasFile, **VoxelizationParameter);
io.writeData(os.path.join(BaseDirectory, 'cells_heatmap_weighted.tif'), vox.astype('int32'));

##################### Label points by anatomical regions

#points = io.readPoints(TransformedCellFile);
counts = countPointsInRegions(points, AnnotationFile);
io.writePoints(os.path.join(BaseDirectory,'Annotated_Points.csv'), counts);

#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################
#####################

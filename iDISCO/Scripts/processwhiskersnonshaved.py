
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 14:58:07 2015

@author: mtllab
"""


import os, numpy, math

import iDISCO.Settings as settings


from iDISCO.Alignment.Resampling import resampleData;

######################### Data parameters

BaseDirectory = '/home/mtllab/Documents/whiskers/nonshaved/150620-6R';
   
cFosFile      = '/home/mtllab/Documents/whiskers/nonshaved/150620-6R/150717_0_8xs3-cfos-20HFcont_06-56-36/06-56-36_0_8xs3-cfos-20HFcont_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
cFosFileRange = {'x' : all, 'y' : all, 'z' : all};

AutofluoFile  = '/home/mtllab/Documents/whiskers/nonshaved/150620-6R/150717_0_8xs3-autofluor_08-33-50/08-33-50_0_8xs3-autofluor_UltraII_C00_xyz-Table Z\d{4}.ome.tif';

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
    "processes" : 4,
   
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
    
ResamplingParameter["source"] = cFosFile;
ResamplingParameter["sink"]   = os.path.join(BaseDirectory, 'cfos_resampled.tif');
    
ResamplingParameter["resolutionSource"] = OriginalResolution;
ResamplingParameter["resolutionSink"]   = ResolutionAffineCFosAutoFluo;

ResamplingParameter["orientation"] = FinalOrientation;
   
resampleData(**ResamplingParameter);


#Files for Auto-fluorescence (acquisition movements correction)
ResamplingParameter["source"] = AutofluoFile;
ResamplingParameter["sink"]   = os.path.join(BaseDirectory, 'autofluo_for_cfos_resampled.tif');

resampleData(**ResamplingParameter);


#Files for Auto-fluorescence (Atlas Registration)
ResamplingParameter["sink"]  =  os.path.join(BaseDirectory, 'autofluo_resampled.tif');
ResamplingParameter["resolutionSink"]  = AtlasResolution;

resampleData(**ResamplingParameter);



### Align cFos and Autofluo

#directory of the alignment result
parameter.Alignment.AlignmentDirectory = os.path.join(basedirectory, 'elastix_cfos_to_auto');
    
#moving and reference images
parameter.Alignment.MovingImage = os.path.join(basedirectory, 'autofluo_for_cfos_resampled.tif');
parameter.Alignment.FixedImage  = os.path.join(basedirectory, 'cfos_resampled.tif');
parameter.Alignment.FixedImageMask = None;


#elastix parameter files for alignment
parameter.Alignment.AffineParameterFile  = os.path.join(PathReg, 'Par0000affine_acquisition.txt');
parameter.Alignment.BSplineParameterFile = None;
runInitializeElastix(parameter)
resultDirectory = runAlignment(parameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory


# Align Autofluo and Atlas

#directory of the alignment result
parameter.Alignment.AlignmentDirectory = os.path.join(basedirectory, 'elastix_auto_to_atlas');
    
#moving and reference images
parameter.Alignment.MovingImage = AtlasFile;
parameter.Alignment.FixedImage  = os.path.join(basedirectory, 'autofluo_resampled.tif');
parameter.Alignment.FixedImageMask = None;

#elastix parameter files for alignment
parameter.Alignment.AffineParameterFile  = os.path.join(PathReg, 'Par0000affine.txt');
parameter.Alignment.BSplineParameterFile = os.path.join(PathReg, 'Par0000bspline.txt');
resultDirectory = runAlignment(parameter);
print "Aligned cfos with autofluo: result directory: %s" % resultDirectory


# Detect Cells


## Cell Detection

# result files for cell coordinates (csv, vtk or ims)
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'cells.csv');
parameter.ImageProcessing.CellIntensityFile = os.path.join(basedirectory, 'intensities.csv');
parameter.DataSource.ImageFile = cFosFile;
runCellDetection(parameter);


### Transform points from Original c-Fos position to autofluorescence

parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'cells.csv')
parameter.ImageProcessing.CellTransformedCoordinateFile = os.path.join(basedirectory, 'cells_transformed_cfos_to_auto.csv');
parameter.Resampling.DataFiles = cFosFile;

#Resolution of the Data (in um / pixel)
parameter.Resampling.ResolutionData = OriginalResolution;

#Resolution of the Reference / Atlas (in um/ pixel)
parameter.Resampling.ResolutionReference = (16, 16, 16);

#Orientation of the Data set compared to Atlas 
parameter.Resampling.Orientation = FinalOrientation;

##Alignment Parameter
#directory of the alignment result
parameter.Alignment.AlignmentDirectory = os.path.join(basedirectory, 'elastix_cfos_to_auto');
    
#moving and reference images
parameter.Alignment.MovingImage = os.path.join(basedirectory, 'autofluo_for_cfos_resampled.tif');
parameter.Alignment.FixedImage  = os.path.join(basedirectory, 'cfos_resampled.tif');
parameter.Alignment.FixedImageMask = None;

runCellCoordinateTransformationToReference(parameter);




### Transform points from autofluorescence to Atlas

##Cells
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'cells_transformed_cfos_to_auto.csv')
parameter.ImageProcessing.CellTransformedCoordinateFile = os.path.join(basedirectory, 'cells_transformed_to_Atlas.csv');

##Resampling Parameter
#Files
parameter.Resampling.DataFiles = os.path.join(basedirectory, 'autofluo_for_cfos_resampled.tif');

#Resolution of the Data (in um / pixel)
parameter.Resampling.ResolutionData = (16, 16, 16);

#Resolution of the Reference / Atlas (in um/ pixel)
parameter.Resampling.ResolutionReference = AtlasResolution;

#Orientation of the Data set compared to Atlas 
parameter.Resampling.Orientation = (1 ,2, 3);

##Alignment Parameter
#directory of the alignment result
parameter.Alignment.AlignmentDirectory = os.path.join(basedirectory, 'elastix_auto_to_atlas');
    
#moving and reference images
parameter.Alignment.MovingImage = AtlasFile;
parameter.Alignment.FixedImage  = os.path.join(basedirectory, 'cfos_resampled.tif');
parameter.Alignment.FixedImageMask = None;

runCellCoordinateTransformationToReference(parameter);


#################### Heat map generation

points = io.readPoints(os.path.join(basedirectory, 'cells_transformed_to_Atlas.csv'));
pv = parameter.Voxelization;
si = dataSize(AtlasFile);
print si;
vox = voxelize(points, si, average = pv.AveragingDiameter, mode = pv.Mode);
io.writeDataStack(os.path.join(basedirectory, 'cells_heatmap.tif'), vox.astype('int32'));


##################### Label points by anatomical regions

pts =  io.readPoints(os.path.join(basedirectory, 'cells_transformed_to_Atlas.csv'));
counts = iDISCO.Analysis.Label.countPointsInRegions(pts, AnnotationFile);
io.writePoints(os.path.join(basedirectory,'Annotated_Points.csv'), counts);

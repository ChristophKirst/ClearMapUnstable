# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 20:23:41 2015

@author: mtllab
"""


# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 14:58:07 2015

@author: mtllab"""


#from parameter1stHaloperidol import *

execfile('parameter1stHaloperidol.py')


import iDISCO.Visualization.Plot as plt

################

#data = io.readData(cFosFile, **cFosFileRange);
#io.writeData(os.path.join(BaseDirectory, 'raw.tif'), data)

dataSource = cFosFile;
x = (400,1000);
y = (400,1000);
z = (500,1500);


points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-new.npy'),  os.path.join(BaseDirectory,  'intensities-new.npy')));
pointsT, intensitiesT = thresholdPoints(points, intensities, threshold = (1000, 65000), row = (1,0));


data = plt.overlayPoints(dataSource, pointsT, x = x, y = y, z = z, pointColor = None);
io.writeData('/home/mtllab/Documents/Haloperidol/1strun/Haloperidol/checkpoints_new.tif', data)


sh = io.pointShiftFromRange(io.dataSize(dataSource), x = x, y = y, z = z)
pointsS = pointsT - sh;


bgFile = '/home/mtllab/Documents/Haloperidol/1strun/Haloperidol/background/background_Z\d{4}.ome.tif';

data2 = plt.overlayPoints(bgFile, pointsS, pointColor = None);
io.writeData('/home/mtllab/Documents/Haloperidol/1strun/Haloperidol/checkpoints_new_label.tif', data2)


hmFile = '/home/mtllab/Documents/Haloperidol/1strun/Haloperidol/hmax/hmax_Z\d{4}.ome.tif';

data3 = plt.overlayPoints(hmFile, pointsS, pointColor = None);
io.writeData('/home/mtllab/Documents/Haloperidol/1strun/Haloperidol/checkpoints_new_hmax.tif', data3)


idsx = numpy.logical_and(pointsS[:,0] > 443, pointsS[:,0] < 447)
idsy = numpy.logical_and(pointsS[:,1] > 83,  pointsS[:,1] < 86)
idsz = numpy.logical_and(pointsS[:,2] > 65,  pointsS[:,2] < 70)

ids = numpy.logical_and(idsx, idsy);
ids = numpy.logical_and(ids,  idsz);

ids.max()


pointsS[ids]
intensitiesT[ids]



intensities = io.readPoints(os.path.join(BaseDirectory, 'intensities-new.npy'))
import matplotlib.pyplot as mpp
mpp.hist(intensities[:,1], bins = 512)



######



############################################################################### 
# RUN
###############################################################################   


### Resample Fluorescent and CFos images
   
#resampleData(**CorrectionResamplingParameterCfos);

#resampleData(**CorrectionResamplingParameterAutoFluo);

#resampleData(**RegistrationResamplingParameter);


### Align cFos and Autofluo

#resultDirectory  = alignData(**CorrectionAlignmentParameter);

### Align Autofluo and Atlas

#resultDirectory  = alignData(**RegistrationAlignmentParameter);


### Detect Cells

## Cell Detection

#detectCells(**ImageProcessingParameter);

## Save spot retained as legitimate cells, based on intensities:

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minthreshold, maxthreshold), row = (1,0));

io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to reference image size

points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);


points = resamplePoints(**CorrectionResamplingPointsParameter);
    
# transform points
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

# resmaple back to original
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);

## Transform points from corrected to registered
# downscale points to reference image size
RegistrationResamplingPointParameter["pointSource"] = points;
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
ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);

table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);

with open(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

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

BaseDirectory = '/home/mtllab/Documents/Haloperidol/1strun/saline';
   
cFosFile = '/home/mtllab/Documents/Haloperidol/1strun/saline/150701_0_8X-cfos_21-19-09/21-19-09_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/Haloperidol/1strun/saline/150701_0_8X-autofluo_23-51-51/23-51-51_0_8X-autofluo_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
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
minthreshold = 20; #remove points detected in the background, based on the filtered intensities
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
AtlasFile      = os.path.join(PathReg, 'half_template_25_right.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right.tif');

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








############################################################################### 
# RUN
###############################################################################   


### Resample Fluorescent and CFos images
   
#resampleData(**CorrectionResamplingParameterCfos);

#resampleData(**CorrectionResamplingParameterAutoFluo);

#resampleData(**RegistrationResamplingParameter);


### Align cFos and Autofluo

#resultDirectory  = alignData(**CorrectionAlignmentParameter);

### Align Autofluo and Atlas

#resultDirectory  = alignData(**RegistrationAlignmentParameter);


### Detect Cells

## Cell Detection

#detectCells(**ImageProcessingParameter);

## Save spot retained as legitimate cells, based on intensities:

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (minthreshold, maxthreshold), row = (1,0));

io.writePoints((os.path.join(BaseDirectory, 'cells.npy'), os.path.join(BaseDirectory,  'intensities.npy')), (points, intensities));


### Transform points from Original c-Fos position to autofluorescence

## Transform points from original to corrected
# downscale points to reference image size

points = io.readPoints(CorrectionResamplingPointsParameter["pointSource"]);


points = resamplePoints(**CorrectionResamplingPointsParameter);
    
# transform points
points = transformPoints(points, transformDirectory = CorrectionAlignmentParameter["resultDirectory"], indices = False, resultDirectory = None);

# resmaple back to original
CorrectionResamplingPointsInverseParameter["pointSource"] = points;
points = resamplePointsInverse(**CorrectionResamplingPointsInverseParameter);

## Transform points from corrected to registered
# downscale points to reference image size
RegistrationResamplingPointParameter["pointSource"] = points;
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
ids, counts = countPointsInRegions(points, labeledImage = AnnotationFile, intensities = intensities, intensityRow = 0);

table = numpy.zeros(ids.shape, dtype=[('id','int64'),('counts','f8'),('name', 'a256')])
table["id"] = ids;
table["counts"] = counts;
table["name"] = labelToName(ids);

with open(os.path.join(BaseDirectory, 'Annotated_counts_intensities.csv'),'w') as f:
    for sublist in table:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();


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

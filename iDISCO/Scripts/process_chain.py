# -*- coding: utf-8 -*-
"""
Main Script to Run Full Pricess Chain

@author: ckirst
"""

# path to illastik installation
sys.path.insert(1, '/usr/local/ilastik-05-rc-final')


from iDISCO.Parameter import * 
import iDISCO.IO.IO as io
import iDISCO.Alignment.Elastix as elx;
import iDISCO.Alignment.Resampling as rspl;
import iDISCO.Analysis.Label as lbl
import iDISCO.Analysis.Statistics as stats

verbose = True;

par = Parameter();


######################### Data parameters

basedirectory = '/home/mtllab/Documents/whiskers/nonshaved/150620-6R';
   
cFosFile = '/home/mtllab/Documents/whiskers/nonshaved/150620-6R/150717_0_8xs3-cfos-20HFcont_06-56-36/06-56-36_0_8xs3-cfos-20HFcont_UltraII_C00_xyz-Table Z\d{4}.ome.tif';
AutofluoFile = '/home/mtllab/Documents/whiskers/nonshaved/150620-6R/150717_0_8xs3-autofluor_08-33-50/08-33-50_0_8xs3-autofluor_UltraII_C00_xyz-Table Z\d{4}.ome.tif';

#Resolution of the Data (in um / pixel)
OriginalResolution = (4.0625, 4.0625, 3);

#Resolution of the Atlas (in um/ pixel)
AtlasResolution = (25, 25, 25);

#Orientation
FinalOrientation = (1,2,3);



############################## Cell Detection Parameters using DoG

parameter.ImageProcessing.Method = "SpotDetection"

#image ranges
parameter.DataSource.XRange = all;
parameter.DataSource.YRange = all;
parameter.DataSource.ZRange = all;

# radius for background removal
parameter.ImageProcessing.Parameter.Background = (15,15);

# size of differeence of gaussian filter
parameter.ImageProcessing.Parameter.Dog = (7, 7, 11);

# h value for h max detection
parameter.ImageProcessing.Parameter.HMax = 10;
parameter.ImageProcessing.Parameter.Threshold = 20;

# intensity threshold for final cells to be saved
parameter.ImageProcessing.Parameter.ThresholdSave = 0;

#################### Heat map generation

##Voxelization

#Size of averaging window
parameter.Voxelization.AveragingDiameter = (15,15,15);  #Radius of the sphere
    
#Image size
parameter.Voxelization.Size = None;  #None extract from data
    
#Mode
parameter.Voxelization.Mode = 'Spherical';

parameter.Voxelization.File = os.path.join(basedirectory, 'points_voxelized.tif');


############################ Config parameters

#Processes to use for Resampling
parameter.Resampling.Processes = 12;

#Elastix binary
parameter.Alignment.ElastixDirectory = '/usr/local/elastix'

#Path to registration parameters and atlases
PathReg = '/home/mtllab/Documents/warping';
AtlasFile = os.path.join(PathReg, 'half_template_25_right.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right.tif');

#Stack Processing Parameter for cell detection using DoG

#max number of parallel processes
parameter.StackProcessing.Processes = 3;
#chunk sizes
parameter.StackProcessing.ChunkSizeMax = 100;
parameter.StackProcessing.ChunkSizeMin = 30;
parameter.StackProcessing.ChunkOverlap = 15;

#optimize chunk size and number to number of processes
parameter.StackProcessing.OptimizeChunks = True;
#increase chunk size for optimizaition (True, False or all = choose automatically)
parameter.StackProcessing.OptimizeChunkSizeIncrease = all;





######################## RUN


## Resample Fluorescent and CFos images
    
parameter.Resampling.ResolutionData = OriginalResolution;
parameter.Resampling.ResolutionReference = (16, 16, 16);
parameter.Resampling.Orientation = FinalOrientation;

#Files for Cfos
parameter.Resampling.DataFiles = cFosFile;
parameter.Resampling.ResampledFile = os.path.join(basedirectory, 'cfos_resampled.tif');
resampledImage = runResampling(parameter);

#Files for Auto-fluorescence (acquisition movements correction)
parameter.Resampling.DataFiles = AutofluoFile;
parameter.Resampling.ResampledFile = os.path.join(basedirectory, 'autofluo_for_cfos_resampled.tif');
resampledImage = runResampling(parameter);


#Files for Auto-fluorescence (Atlas Registration)
parameter.Resampling.ResampledFile = os.path.join(basedirectory, 'autofluo_resampled.tif');
parameter.Resampling.ResolutionReference = AtlasResolution;
resampledImage = runResampling(parameter);


# Align cFos and Autofluo

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

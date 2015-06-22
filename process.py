# -*- coding: utf-8 -*-
"""
Processing Brain Samples - Main Example Script


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import os

from iDISCO.Parameter import *
from iDISCO.Run import *
from iDISCO.IO import IO as io


##################################
# Segmentation
##################################

parameter = Parameter();

basedirectory = iDISCOPath()

### Data Source 

#raw data from microscope used for cell detection (ims or tif)
parameter.DataSource.ImageFile = os.path.join(basedirectory, 'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z\d{4}.ome.tif')

#image ranges
parameter.DataSource.XRange = all;
parameter.DataSource.YRange = all;
parameter.DataSource.ZRange = all;


### Cell Segmentation Parameter

# method to segment cells (SpotDetection or Ilastik)
parameter.ImageProcessing.Method = 'SpotDetection';

# radius for background removal
parameter.ImageProcessing.Parameter.Background = (15,15);

# size of differeence of gaussian filter
parameter.ImageProcessing.Parameter.Dog = (11, 7, 7);

# h value for h max detection
parameter.ImageProcessing.Parameter.HMax = 20;

# intensity threshold for cells 
parameter.ImageProcessing.Parameter.Threshold = 20;

# intensity threshold for final cells to be saved
parameter.ImageProcessing.Parameter.ThresholdSave = 30;

# result file for cell coordinates (csv, vtk or ims)
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'Test/ImageProcessing/cells.csv');


### Parallel Processing Parameter
    
#max number of parallel processes
parameter.ParallelProcessing.Processes = 2;
   
#chunk sizes
parameter.ParallelProcessing.ChunkSizeMax = 100;
parameter.ParallelProcessing.ChunkSizeMin = 30;
parameter.ParallelProcessing.ChunkOverlap = 15;

#optimize chunk size and number to number of processes
parameter.ParallelProcessing.OptimizeChunks = True;
    
#increase chunk size for optimizaition (True, False or all = choose automatically)
parameter.ParallelProcessing.OptimizeChunkSizeIncrease = all;



### Run and Save

centers, intensities = runCellDetection(parameter);

iid = intensities > parameter.ImageProcessing.Parameter.ThresholdSave ;
centersSave = centers[iid,:];

io.writePoints(parameter.ImageProcessing.CellCoordinateFile, centersSave);



##################################
# Transform Cell Coordinates
##################################

parameter = Parameter();

basedirectory = iDISCOPath()

### Alignment Parameter

#Transformix binary
parameter.Alignment.Transformix = '/home/ckirst/programs/elastix/bin/transformix'

#Directory of elastix output form Alignment
parameter.Alignment.AlignmentDirectory = os.path.join(basedirectry, 'Test/Elastix/Output')

#Cell Coordinates csv file
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'Test/ImageProcessing/cells.csv');

#File for transfomred coordinates (csv, vtk or ims)
parameter.Alignment.CellCoordinateFile =  os.path.join(basedirectory, 'Test/Elastix/cells.csv');

# run
points = runCellCoordinateTransformation(parameter);

#save 
io.writePoints(parameter.Alignment.CellCoordinateFile, points);





##################################
# Align Data
##################################

parameter = Parameter();

basedirectory = iDISCOPath()

### Alignment Parameter

#Elastix binary
parameter.Alignment.Elastix = '/home/ckirst/programs/elastix/bin/elastix'

#Directory of elastix output form Alignment
parameter.Alignment.AlignmentDirectory = os.path.join(basedirectry, '/Test/Elastix/Output')

#moving and reference images
parameter.Alignment.MovingImage = os.path.join(basedirectry, 'Test/Data/Elastix/150524_0_8X-s3-20HFautofluor_18-51-1-warpable.tif');
parameter.Alignment.FixedImage  = os.path.join(basedirectry, 'Test/Data/Elastix/OstenRefARA_v2_lowerHalf.tif');

#elastix parameter files for alignment
parameter.Alignment.AffineParameterFile  = os.path.join(basedirectry, 'Test/Elastix/ElastixParamterAffine.txt');
parameter.Alignment.BSplineParameterFile = os.path.join(basedirectry, 'Test/Elastix/ElastixParamterBSpline.txt');

#run alignment
runAlignment(parameter);








##################################
# Segmentation using Ilastik
##################################

parameter = Parameter();

basedirectory = iDISCOPath()

### Data Source 

#raw data from microscope used for cell detection (ims or tif)
parameter.DataSource.ImageFile = os.path.join(basedirectory, 'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z\d{4}.ome.tif')

#image ranges
parameter.DataSource.XRange = all;
parameter.DataSource.YRange = all;
parameter.DataSource.ZRange = all;


### Cell Segmentation Parameter

# method to segment cells (SpotDetection or Ilastik)
parameter.ImageProcessing.Method = 'Ilastik';

# radius for background removal
parameter.ImageProcessing.Parameter.Background = (15,15);

# result file for cell coordinates (csv, vtk or ims)
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'Test/ImageProcessing/cells.csv');


### Run and Save

centers, intensities = runCellDetection(parameter);

iid = intensities > parameter.ImageProcessing.Parameter.ThresholdSave ;
centersSave = centers[iid,:];

io.writePoints(parameter.ImageProcessing.CellCoordinateFile, centersSave);















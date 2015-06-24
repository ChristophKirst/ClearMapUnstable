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
parameter.DataSource.ImageFile = '/run/media/ckirst/ChristophsBackuk4TB/Data/Science/Projects/BrainActiityMap/Experiment/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims'
#image ranges
parameter.DataSource.XRange = (1400, 1800);
parameter.DataSource.YRange = (1000, 1300);
parameter.DataSource.ZRange = (300, 600);


### Cell Segmentation Parameter

# method to segment cells (SpotDetection or Ilastik)
parameter.ImageProcessing.Method = 'SpotDetection';

# radius for background removal
parameter.ImageProcessing.Parameter.Background = (15,15);

# size of differeence of gaussian filter
parameter.ImageProcessing.Parameter.Dog = (7, 7, 11);

# h value for h max detection
parameter.ImageProcessing.Parameter.HMax = 20;

# intensity threshold for cells 
parameter.ImageProcessing.Parameter.Threshold = 20;

# intensity threshold for final cells to be saved
parameter.ImageProcessing.Parameter.ThresholdSave = 30;

# result file for cell coordinates (csv, vtk or ims)
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'Test/ImageProcessing/cells.csv');


### Stack Processing Parameter
    
#max number of parallel processes
parameter.StackProcessing.Processes = 2;
   
#chunk sizes
parameter.StackProcessing.ChunkSizeMax = 100;
parameter.StackProcessing.ChunkSizeMin = 30;
parameter.StackProcessing.ChunkOverlap = 15;

#optimize chunk size and number to number of processes
parameter.StackProcessing.OptimizeChunks = True;
    
#increase chunk size for optimizaition (True, False or all = choose automatically)
parameter.StackProcessing.OptimizeChunkSizeIncrease = all;



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

parameter.Alignment.Lib = '/home/ckirst/programs/elastix/lib'

#Directory of elastix output form Alignment
parameter.Alignment.AlignmentDirectory = os.path.join(basedirectory, 'Test/Elastix/Output')

#Cell Coordinates csv file
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'Test/ImageProcessing/cells.csv');

#File for transfomred coordinates (csv, vtk or ims)
parameter.Alignment.CellCoordinateFile =  os.path.join(basedirectory, 'Test/Elastix/cells_aligned.csv');

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

parameter.Alignment.Lib = '/home/ckirst/programs/elastix/lib'

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

import os

from iDISCO.Parameter import *
from iDISCO.Run import *
from iDISCO.IO import IO as io
import math


parameter = Parameter();

basedirectory = iDISCOPath()

### Data Source 

#raw data from microscope used for cell detection (ims or tif)
parameter.DataSource.ImageFile = os.path.join(basedirectory, 'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z\d{4}.ome.tif')

#image ranges
parameter.DataSource.XRange = (1400, 1800);
parameter.DataSource.YRange = (1000, 1300);
parameter.DataSource.ZRange = all;


### Cell Segmentation Parameter

# method to segment cells (SpotDetection or Ilastik)
parameter.ImageProcessing.Method = 'Ilastik';

parameter.ImageProcessing.Parameter = IlastikParameter()

#rescale to fi uint8 data type and enhanve contrast for ilastik
parameter.ImageProcessing.Parameter.Rescale  =  1.0 / math.pow(2,16) * math.pow(2,8) * 10;

# radius for background removal
parameter.ImageProcessing.Parameter.Background = (15,15);

# classifier
parameter.ImageProcessing.Parameter.Classifier = os.path.join(basedirectory, "Test/Ilastik/classifier.h5");

# result file for cell coordinates (csv, vtk or ims)
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'Test/ImageProcessing/cells.csv');


#Stack Processing
  
#chunk sizes
parameter.StackProcessing.ChunkSizeMax = 100;
parameter.StackProcessing.ChunkSizeMin = 30;
parameter.StackProcessing.ChunkOverlap = 15;


### Run and Save

centers, intensities = runCellDetection(parameter);

iid = intensities > parameter.ImageProcessing.Parameter.ThresholdSave ;
centersSave = centers[iid,:];

io.writePoints(parameter.ImageProcessing.CellCoordinateFile, centersSave);















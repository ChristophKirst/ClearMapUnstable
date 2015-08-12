# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 20:27:23 2015

@author: mtllab
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 14:58:07 2015

@author: mtllab
"""


import os
import numpy
from iDISCO.Parameter import * 
from iDISCO.Run import *
from iDISCO.Analysis.Label import *

verbose = True;

parameter = Parameter();

######################### Data parameters

basedirectory = '/home/mtllab/Documents/countingaccuracy';
   
cFosFile = '/home/mtllab/Documents/countingaccuracy/countingaccuracy.tif';


#Resolution of the Data (in um / pixel)
OriginalResolution = (4.0625, 4.0625, 3);

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
parameter.ImageProcessing.Parameter.HMax = 100;
parameter.ImageProcessing.Parameter.Threshold = 20;

# intensity threshold for final cells to be saved
parameter.ImageProcessing.Parameter.ThresholdSave = 0;

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



## Cell Detection

# result files for cell coordinates (csv, vtk or ims)
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'cells.csv');
parameter.ImageProcessing.CellIntensityFile = os.path.join(basedirectory, 'intensities.csv');
parameter.DataSource.ImageFile = cFosFile;
runCellDetection(parameter);

#Point writing as Tiff

from iDISCO.Alignment.Resampling import dataSize
import iDISCO.Analysis.Voxelization as vox

pts = io.readPoints(parameter.ImageProcessing.CellCoordinateFile);
ds = dataSize(cFosFile);
ds = (ds[1], ds[0], ds[2]);

voximg = vox.voxelizePixel(pts, ds);
io.writeData(os.path.join(basedirectory, 'points_pixel.tif'), voximg)


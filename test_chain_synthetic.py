# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 11:42:52 2015

Test Processing Chain with Synthetic Data

@author: ckirst
"""



##############################################################################
# Test Resampling / Make Reference
##############################################################################

import os
import numpy

from iDISCO.Parameter import *
from iDISCO.Run import runResampling

import shutil

basedirectory = os.path.join(iDISCOPath(), 'Test');

parameter = Parameter();


#Files
parameter.Resampling.DataFiles = os.path.join(basedirectory, 'Data/Synthetic/test_iDISCO_\d{3}.tif')
parameter.Resampling.ResampledFile = os.path.join(basedirectory, 'Synthetic/test_iDISCO_resample.tif');
    
#Resolution of the Data (in um / pixel)
parameter.Resampling.ResolutionData = (5, 5, 3);

#Resolution of the Reference / Atlas (in um/ pixel)
parameter.Resampling.ResolutionReference = (12, 15, 5);

#Orientation of the Data set wrt reference 
#(-axis will invert the orientation, for other hemisphere use (-1, 2, 3), to exchnge x,y use (2,1,3) etc)
parameter.Resampling.Orientation = (1,2,3);

#Processes to use for Resmapling
parameter.Resampling.Processes = 4;

resampledImage = runResampling(parameter);

print "Resampled image saved as %s" % resampledImage

# Reference = Identity to test
shutil.copyfile(resampledImage, os.path.join(basedirectory, 'Synthetic/test_iDISCO_reference.tif'))



##############################################################################
# Test Alignment
##############################################################################

import os

from iDISCO.Parameter import *
from iDISCO.Run import runAlignment, runInitializeElastix

import iDISCO.Alignment.Elastix

basedirectory = os.path.join(iDISCOPath(), 'Test');

parameter = Parameter();

#directory of the alignment result
parameter.Alignment.AlignmentDirectory = os.path.join(basedirectory, 'Synthetic/elastix');

#Elastix binary
parameter.Alignment.ElastixDirectory = '/home/ckirst/programs/elastix'
    
#moving and reference images
parameter.Alignment.MovingImage = os.path.join(basedirectory, 'Synthetic/test_iDISCO_resample.tif');
parameter.Alignment.FixedImage  = os.path.join(basedirectory, 'Synthetic/test_iDISCO_reference.tif');
parameter.Alignment.FixedImageMask = None;
  
  
#elastix parameter files for alignment
parameter.Alignment.AffineParameterFile  = os.path.join(parameter.Alignment.AlignmentDirectory, 'ElastixParameterAffine.txt');
#parameter.Alignment.BSplineParameterFile = os.path.join(parameter.Alignment.AlignmentDirectory, 'ElastixParameterBSpline.txt');
parameter.Alignment.BSplineParameterFile = None;


runInitializeElastix(parameter)
iDISCO.Alignment.Elastix.ElastixSettings.printInfo();

resultDirectory = runAlignment(parameter);

print "Aligned images: result directory: %s" % resultDirectory




##############################################################################
# Test Detect Points
############################################################################## 

import os

from iDISCO.Parameter import *
from iDISCO.Run import runCellDetection

parameter = Parameter();

basedirectory = os.path.join(iDISCOPath(), 'Test');

### Data Source 

#raw data from microscope used for cell detection (ims or tif)
parameter.DataSource.ImageFile = os.path.join(basedirectory, 'Data/Synthetic/test_iDISCO_\d{3}.tif')

#image ranges
parameter.DataSource.XRange = all;
parameter.DataSource.YRange = all;
parameter.DataSource.ZRange = all;

# radius for background removal
parameter.ImageProcessing.Parameter.Background = (20,20);

# size of differeence of gaussian filter
parameter.ImageProcessing.Parameter.Dog = (8, 8, 6);

# h value for h max detection
parameter.ImageProcessing.Parameter.HMax = 0.5;
parameter.ImageProcessing.Parameter.Threshold = 0.5;

# intensity threshold for final cells to be saved
parameter.ImageProcessing.Parameter.ThresholdSave = None;

# result file for cell coordinates (csv, vtk or ims)
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'Synthetic/cells.csv');
parameter.ImageProcessing.CellIntensityFile = os.path.join(basedirectory, 'Synthetic/intensities.csv');

runCellDetection(parameter)



##############################################################################
# Test Resample Points
############################################################################## 

import os
import numpy

from iDISCO.Parameter import *
from iDISCO.Run import runResampling
import iDISCO.IO.IO as io
from iDISCO.Visualization import Plot

from iDISCO.Alignment.Resampling import resamplePoints, resamplePointsInverse, dataSize

import shutil

verbose = True;

basedirectory = os.path.join(iDISCOPath(), 'Test');

parameter = Parameter();

files = os.path.join(basedirectory, 'Data/Synthetic/test_iDISCO_\d{3}.tif');

datasize = dataSize(files);
print datasize

#Resolution of the Data (in um / pixel)
parameter.Resampling.ResolutionData = (5, 5, 3);

#Resolution of the Reference / Atlas (in um/ pixel)
parameter.Resampling.ResolutionReference = (12, 15, 5);

#Orientation of the Data set wrt reference 
#(-axis will invert the orientation, for other hemisphere use (-1, 2, 3), to exchnge x,y use (2,1,3) etc)
parameter.Resampling.Orientation = (1,2,3);

#centers
centers = io.readPoints(os.path.join(basedirectory, 'Synthetic/cells.csv'))
print "Center shape: (%d, %d)" % centers.shape
rcenters = resamplePoints(centers, datasize, parameter.Resampling.ResolutionData, parameter.Resampling.ResolutionReference, parameter.Resampling.Orientation)
print "Reshaped centers shape: (%d, %d)" % rcenters.shape

if verbose:
    dataraw = io.readData(os.path.join(basedirectory, 'Data/Synthetic/test_iDISCO_\d{3}.tif'));
    datares = io.readData(os.path.join(basedirectory, 'Synthetic/test_iDISCO_resample.tif'));
    print "Shape raw: " + str(dataraw.shape)
    print "Shape res: " + str(datares.shape)
    
    Plot.plotOverlayPoints(dataraw*0.01, centers)
    Plot.plotOverlayPoints(datares*0.01, rcenters)


#check inverse:
rrcenters = resamplePointsInverse(rcenters, datasize, parameter.Resampling.ResolutionData, parameter.Resampling.ResolutionReference, parameter.Resampling.Orientation);
diff = centers - rrcenters;
print (diff.max(), diff.min());



##############################################################################
# Test Transform Points
############################################################################## 

import iDISCO.Visualization.Plot as Plot;

import iDISCO.Alignment.Elastix as elx;

elx.initializeElastix('/home/ckirst/programs/elastix')
elx.ElastixSettings.printInfo()

verbose = True;

transformdir = os.path.join(basedirectory, 'Synthetic/transformix')

pointsdir = os.path.join(basedirectory, 'Synthetic/points')

parameter.Alignment.AlignmentDirectory = os.path.join(basedirectory, 'Synthetic/elastix');

acenters = elx.transformPoints(rcenters, alignmentdirectory = parameter.Alignment.AlignmentDirectory, outdirectory = pointsdir);
acenters = elx.transformPoints(rcenters, alignmentdirectory = parameter.Alignment.AlignmentDirectory);


if verbose:
    transformdata = io.readData(os.path.join(transformdir, 'result.mhd'));
    referencedata = io.readData(os.path.join(basedirectory, 'Synthetic/test_iDISCO_resample.tif'))
    Plot.plotOverlayPoints(transformdata * 0.01, acenters);
    Plot.plotOverlayPoints(referencedata * 0.01, acenters);


##############################################################################
# Test Transform Resampled Data to Reference
############################################################################## 

from iDISCO.Parameter import *
import iDISCO.Alignment.Elastix as elx
import iDISCO.Visualization.Plot as Plot
import iDISCO.IO.IO as io

basedirectory = os.path.join(iDISCOPath(), 'Test');

verbose = True;

resultdir = os.path.join(basedirectory, 'Synthetic/elastix');
transformfile = elx.getTransformParameterFile(resultdir)

# tranform points according to alignment of data and reference 
resamplefile = os.path.join(basedirectory, 'Synthetic/test_iDISCO_resample.tif');
transformdir = os.path.join(basedirectory, 'Synthetic/transformix')

elx.initializeElastix('/home/ckirst/programs/elastix')
elx.ElastixSettings.printInfo()

#elx.transformData(dataresname, alignmentdirectory = resultdir, outdirectory = dataalgname)
elx.transformData(resamplefile, transformparameterfile = os.path.join(basedirectory, 'Synthetic/elastix/TransformParameters.0.txt'), outdirectory = transformdir)

if verbose:
    transformdata = io.readData(os.path.join(transformdir, 'result.mhd'));
    resampledata = io.readData(resamplefile);

    Plot.plotTiling(0.01 * transformdata)
    Plot.plotTiling(0.01 * resampledata)
    










# transform raw points to aligned points in raw data size









# voxelize and plot




##############################################################################
# Test Full Chain
##############################################################################

parameter = Parameter();

basedirectory = iDISCOPath()










### Data Source 

#raw data from microscope used for cell detection (ims or tif)
parameter.DataSource.ImageFile = os.path.join(basedirectory, 'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z\d{4}.ome.tif')

parameter.DataSource.ImageFile = '/run/media/ckirst/ChristophsBackuk4TB/Data/Science/Projects/BrainActiityMap/Experiment/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims'

#image ranges
parameter.DataSource.XRange = all;
parameter.DataSource.YRange = all;
parameter.DataSource.ZRange = (100,120);

# load data

data = io.readData(parameter.DataSource.ImageFile, x = parameter.DataSource.XRange, y = parameter.DataSource.YRange, z = parameter.DataSource.ZRange, resolution = 0);
print "Loaded data from " + parameter.DataSource.ImageFile;
print "Data size is: " + str(data.shape)

# visualize

if verbose:
    plt.plotTiling(15*data)


### Process using Spot Detection

# radius for background removal
parameter.ImageProcessing.Parameter.Background = (15,15);

img = ip.removeBackground(data, parameter = parameter.ImageProcessing, verbose =  verbose);
print img.dtype
print img.shape

# size of differeence of gaussian filter
parameter.ImageProcessing.Parameter.Dog = (7, 7, 11);

img = ip.dogFilter(img, parameter = parameter.ImageProcessing, verbose =  verbose);
print img.dtype
print img.shape

# h value for h max detection
parameter.ImageProcessing.Parameter.HMax = 20;

imgmax = ip.findExtendedMaxima(img, parameter = parameter.ImageProcessing, verbose =  verbose);
print img.dtype
print img.shape

centers, intensities = ip.findCenterOfMaxima(img, imgmax, parameter = parameter.ImageProcessing, verbose =  verbose);


# intensity threshold for final cells to be saved
parameter.ImageProcessing.Parameter.ThresholdSave = 30;

iid = intensities > parameter.ImageProcessing.Parameter.ThresholdSave ;
centersSave = centers[iid,:];

if verbose:
    imgc = numpy.zeros(img.shape);
    for i in range(centersSave.shape[0]):
        imgc[centersSave[i,0], centersSave[i,1], centersSave[i,2]] = 1;
    plt.plotOverlayLabel(img * 0.01, imgc, alpha = False);



# result file for cell coordinates (csv, vtk or ims)
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'Test/ImageProcessing/cells.csv');

io.writePoints(parameter.ImageProcessing.CellCoordinateFile, centersSave);






















##############################################################################
# Test Spot Detection
##############################################################################

import os
import numpy
import sys

from iDISCO.Parameter import *
from iDISCO.IO import IO as io
from iDISCO.Visualization import Plot as plt
from iDISCO.ImageProcessing import SpotDetection as ip

verbose = True;


parameter = Parameter();

basedirectory = iDISCOPath()
print basedirectory

### Data Source 

#raw data from microscope used for cell detection (ims or tif)
parameter.DataSource.ImageFile = os.path.join(basedirectory, 'Test/Data/Synthetic/test_iDISCO_\d{3}.tif')

#image ranges
parameter.DataSource.XRange = all;
parameter.DataSource.YRange = all;
parameter.DataSource.ZRange = all;

# load data

data = io.readData(parameter.DataSource.ImageFile, x = parameter.DataSource.XRange, y = parameter.DataSource.YRange, z = parameter.DataSource.ZRange, resolution = 0);
print "Loaded data from " + parameter.DataSource.ImageFile;
print "Data size is: " + str(data.shape)

# visualize

if verbose:
    plt.plotTiling(15*data)

### Process using Spot Detection

# radius for background removal
parameter.ImageProcessing.Parameter.Background = (15,15);

img = ip.removeBackground(data, parameter = parameter.ImageProcessing, verbose = verbose);
print img.dtype
print img.shape

# size of differeence of gaussian filter
parameter.ImageProcessing.Parameter.Dog = (7, 7, 5);

img = ip.dogFilter(img, parameter = parameter.ImageProcessing, verbose =  verbose);
print img.dtype
print img.shape

plt.plotTiling(img)


# h value for h max detection
parameter.ImageProcessing.Parameter.HMax = 1;
parameter.ImageProcessing.Parameter.Threshold = 1;

imgmax = ip.findExtendedMaxima(img, parameter = parameter.ImageProcessing, verbose =  False);
print imgmax.max()

centers, intensities = ip.findCenterOfMaxima(img, imgmax, verbose =  verbose);


# intensity threshold for final cells to be saved
parameter.ImageProcessing.Parameter.ThresholdSave = 1;

iid = intensities > parameter.ImageProcessing.Parameter.ThresholdSave ;
centersSave = centers[iid,:];

if verbose:
    imgc = numpy.zeros(img.shape);
    for i in range(centersSave.shape[0]):
        imgc[centersSave[i,0], centersSave[i,1], centersSave[i,2]] = 1;
    plt.plotOverlayLabel(data[:,:,15:30] * 0.1, imgc[:,:,15:30], alpha = False);



# result file for cell coordinates (csv, vtk or ims)
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'Test/ImageProcessing/cells.csv');

io.writePoints(parameter.ImageProcessing.CellCoordinateFile, centersSave);


















# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 16:13:56 2015

@author: ckirst
"""


##############################################################################
# Test Alignment
##############################################################################

import os

from iDISCO.Parameter import *
from iDISCO.Run import runAlignment, runInitializeElastix

import iDISCO.Alignment.Elastix

basedirectory = '/home/ckirst/Science/Projects/BrainActivityMap/Experiment/FullChainTest/';

parameter = Parameter();

#directory of the alignment result
parameter.Alignment.AlignmentDirectory = os.path.join(basedirectory, 'elastix');

#Elastix binary
parameter.Alignment.ElastixDirectory = '/home/ckirst/programs/elastix'
    
#moving and reference images
parameter.Alignment.MovingImage = os.path.join(basedirectory, 'half_template_25_right.tif');
parameter.Alignment.FixedImage  = os.path.join(basedirectory, 'autofluo_resample.tif');
parameter.Alignment.FixedImageMask = None;
  
  
#elastix parameter files for alignment
pp = '/home/ckirst/Science/Projects/BrainActivityMap/Analysis/iDISCO/Test/Elastix/'
parameter.Alignment.AffineParameterFile  = os.path.join(pp, 'ElastixParameterAffine.txt');
parameter.Alignment.BSplineParameterFile = os.path.join(pp, 'ElastixParameterBSpline.txt');
#parameter.Alignment.BSplineParameterFile = None;


runInitializeElastix(parameter)
iDISCO.Alignment.Elastix.ElastixSettings.printInfo();

resultDirectory = runAlignment(parameter);

print "Aligned images: result directory: %s" % resultDirectory




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


verbose = False;

basedirectory = '/home/ckirst/Science/Projects/BrainActivityMap/Experiment/FullChainTest/';

parameter = Parameter();

#Resolution of the Data (in um / pixel)
parameter.Resampling.ResolutionData = (4.0625, 4.0625, 3);

#Resolution of the Reference / Atlas (in um/ pixel)
parameter.Resampling.ResolutionReference = (25, 25, 25);

#Orientation of the Data set wrt reference
#(-axis will invert the orientation, for other hemisphere use (-1, 2, 3), to exchnge x,y use (2,1,3) etc)
parameter.Resampling.Orientation = (1,2,3);

#centers
centers = io.readPoints(os.path.join(basedirectory, 'cells.csv'))
print "Center shape: (%d, %d)" % centers.shape


datasize = (2560, 2160, 1728); # takes (y,x,z)
print datasize

rcenters = resamplePoints(centers, datasize, parameter.Resampling.ResolutionData, parameter.Resampling.ResolutionReference, parameter.Resampling.Orientation)
print "Reshaped centers shape: (%d, %d)" % rcenters.shape

if verbose:
    datares = io.readData(os.path.join(basedirectory, 'autofluo_resample.tif'));
    #print "Shape raw: " + str(dataraw.shape)
    print "Shape res: " + str(datares.shape)
    
    #Plot.plotOverlayPoints(dataraw*0.01, centers)
    Plot.plotOverlayPoints(datares*0.01, rcenters)
    
io.writePoints(os.path.join(basedirectory, 'cells_resampled.csv'), rcenters)

#save pixel image:
import iDISCO.Analysis.Voxelization as vox;
import iDISCO.IO.IO as io


rdata = io.readData(os.path.join(basedirectory, 'autofluo_resample.tif'));
rcenters = io.readPoints(os.path.join(basedirectory, 'cells_resampled.csv'));
voximg = vox.voxelizePixel(rcenters, rdata.shape);
io.writeDataStack(os.path.join(basedirectory, 'cells_resampled.tif'), voximg.astype('int16'))




##############################################################################
# Transform Points from Raw Data to Reference
############################################################################## 

import os

from iDISCO.Parameter import *
from iDISCO.Run import runInitializeElastix, runCellCoordinateTransformationToReference


import iDISCO.Visualization.Plot as Plot
import iDISCO.IO.IO as io

import iDISCO.Analysis.Voxelization as vox;


basedirectory = '/home/ckirst/Science/Projects/BrainActivityMap/Experiment/FullChainTest/';

verbose = True;

parameter = Parameter();

##Cells
parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'cells.csv')
parameter.ImageProcessing.CellTransformedCoordinateFile = os.path.join(basedirectory, 'cells_transformed_1.csv');
##Resampling Parameter
#Files
parameter.Resampling.DataFiles = (2160, 2560, 1728);
parameter.Resampling.ResampledFile = os.path.join(basedirectory, 'Synthetic/test_iDISCO_resample.tif');


#Resolution of the Data (in um / pixel)
parameter.Resampling.ResolutionData = (4.0625, 4.0625, 3);

#Resolution of the Reference / Atlas (in um/ pixel)
parameter.Resampling.ResolutionReference = (25, 25, 25);

#Orientation of the Data set wrt reference
#(-axis will invert the orientation, for other hemisphere use (-1, 2, 3), to exchnge x,y use (2,1,3) etc)
parameter.Resampling.Orientation = (1,2,3);


##Alignment Parameter
#directory of the alignment result
parameter.Alignment.AlignmentDirectory = os.path.join(basedirectory, 'elastix');

#Elastix binary
parameter.Alignment.ElastixDirectory = '/home/ckirst/programs/elastix'
    
#moving and reference images
parameter.Alignment.MovingImage = os.path.join(basedirectory, 'autofluo_resample.tif');
parameter.Alignment.FixedImage  = os.path.join(basedirectory, 'half_template_25_right.tif');
parameter.Alignment.FixedImageMask = None;
  
#elastix parameter files for alignment
#parameter.Alignment.AffineParameterFile  = os.path.join(parameter.Alignment.AlignmentDirectory, '');
#parameter.Alignment.BSplineParameterFile = os.path.join(parameter.Alignment.AlignmentDirectory, 'ElastixParameterBSpline.txt');#
#parameter.Alignment.BSplineParameterFile = None;

runInitializeElastix(parameter)


parameter.ImageProcessing.CellCoordinateFile = os.path.join(basedirectory, 'cells.csv')
pts = io.readPoints(parameter.ImageProcessing.CellCoordinateFile);
#pts = pts[0:10000, :];
#pts[0,:] = [0,0,0];
#pts = pts[:, [1,0,2]];

parameter.ImageProcessing.CellCoordinateFile = pts;

runCellCoordinateTransformationToReference(parameter)

pts2 = io.readPoints(parameter.ImageProcessing.CellTransformedCoordinateFile);
#pts2 = pts2[:,[1,0,2]];

if verbose:
    refdata = io.readData(parameter.Alignment.FixedImage);
    #pts = io.readPoints(parameter.ImageProcessing.CellTransformedCoordinateFile);
    Plot.plotOverlayPoints(0.1 * refdata, pts2)




pts2 = io.readPoints(parameter.ImageProcessing.CellTransformedCoordinateFile);
#pts2 = pts2[:,[1,0,2]];
voximg = vox.voxelizePixel(pts2, refdata.shape);
io.writeDataStack(os.path.join(basedirectory, 'points_transformed_pixel.tif'), voximg)


















##############################################################################
# Misc / Test Transform Points
############################################################################## 

import iDISCO.Visualization.Plot as Plot;

import iDISCO.Alignment.Elastix as elx;

import iDISCO.Analysis.Voxelization as vox
import iDISCO.IO.IO as io

elx.initializeElastix('/home/ckirst/programs/elastix')
elx.ElastixSettings.printInfo()

verbose = True;

parameter.Alignment.AlignmentDirectory = os.path.join(basedirectory, 'elastix');

#rcenters = numpy.array([[0,0,0],[0,52,0], [0,0,23], [74,0,0],[74,52,0]]);

#acenters = elx.transformPoints(rcenters, alignmentdirectory = parameter.Alignment.AlignmentDirectory, outdirectory = pointsdir);

#rcenters = io.readPoints(os.path.join(basedirectory, 'cells_resampled.csv'));

rcenters2 = rcenters.copy();
rcenters2 = rcenters2[:, [1,0,2]];
print rcenters2.shape

acenters = elx.transformPoints(rcenters2, alignmentdirectory = parameter.Alignment.AlignmentDirectory);

acenters2 = acenters.copy();
acenters2 = acenters2[:, [1,0,2]];

#print acenters

#rdata = io.readData(os.path.join(basedirectory, 'autofluo_resample.tif'));
rdata = io.readData(os.path.join(basedirectory, 'half_template_25_right.tif'));

voximg = vox.voxelizePixel(acenters2, rdata.shape) * 5000;
io.writeDataStack(os.path.join(basedirectory, 'cells_transformed.tif'), voximg.astype('int16'))



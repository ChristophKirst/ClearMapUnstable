# -*- coding: utf-8 -*-

"""
Example Parameter file

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import os
import sys
self = sys.modules[__name__];

def iDISCOPath():
    """Returns path of IDISCO software"""
    fn = os.path.split(self.__file__)
    fn = os.path.abspath(fn[0]);
    return fn;




##############################################################################
# Image Processing
##############################################################################

class DataSourceParameter(object):
    """Parameter for files"""
    
    #File name or file pattern of raw image data
    ImageFile = None;
         
    #Co-ordinate ranges of raw data
    ZRange = all;
    XRange = all;
    YRange = all;


class SpotDetectionParameter(object):
    """Paramters used by the spot detection algorythm"""
    
    # Background correctoin: None or (x,y) which is size of disk for gray scale opening
    Background = (15,15);
    
    # Spot Detection via Difference of Gaussians (DoG) filter: (x,y,z) size
    Dog = (11, 7, 7);
    
    #h of h-max transform
    HMax = 20;
    
    #Threshold for min intensity at center to be counted as cell (should be similar to the h max)
    Threshold = 20;
    
    
class IlastikParameter(object):
    """Paramters for Ilastik classification"""
    
    #ilastik 0.5 path
    Ilastik = None;
    
    #ilastic classifier to use
    Classifier = os.path.join(self.iDISCOPath(), '/Test/Ilastik/classifier.h5');  
    
    # Rescaling of images
    Rescale = None;
    
    # Background correctoin: None or (x,y) which is size of disk for gray scale opening
    Background = (15,15);


class ImageProcessingParameter(object):
    """Parameter for image processing chain"""
    
    Method = 'SpotDetection';
    
    Parameter = SpotDetectionParameter();
    
    #File name for cell coordinates csv, vtk or  ims extension
    CellCoordinateFile = None;
    CellInstensityFile = None;


class StackProcessingParameter(object):
    """Parameter for processing a stack in parallel"""
    
    #max number of parallel processes
    Processes = 2;
   
    #chunk sizes
    ChunkSizeMax = 100;
    ChunkSizeMin = 30;
    ChunkOverlap = 15;

    #optimize chunk size and number to number of processes
    OptimizeChunks = True;
    
    #increase chunk size for optimizaition (True, False or all = choose automatically)
    OptimizeChunkSizeIncrease = all;


##############################################################################
# Alignment / Morphin / Resmapling
##############################################################################

class AlignmentParameter(object):
    """Parameter for Elastix alignment"""
    
    #directory of the alignment result
    AlignmentDirectory = None;
    
    #Elastix binary
    ElastixDirectory = '/usr/local/elastix/'
        
    #moving and reference images
    MovingImage = os.path.join(self.iDISCOPath(), '/Test/Data/Elastix/150524_0_8X-s3-20HFautofluor_18-51-1-warpable.tif');
    FixedImage  = os.path.join(self.iDISCOPath(), '/Test/Data/Elastix/OstenRefARA_v2_lowerHalf.tif');
    FixedImageMask = None;
    
    #elastix parameter files for alignment
    AffineParameterFile  = os.path.join(self.iDISCOPath(), '/Test/Elastix/ElastixParameterAffine.txt');
    BSplineParameterFile = os.path.join(self.iDISCOPath(), '/Test/Elastix/ElastixParameterBSpline.txt');
    

class ResamplingParameter(object):
    """Parameter for resampling data to reference atlas"""
    
    #Data ource and output file
    DataFiles = None;
    ResampledFile = None;
    
    #Resolution of the raw data (in um / pixel)
    ResolutionData = (4.0625, 4.0625, 3);

    #Resolution of the Reference / Atlas (in um/ pixel)
    ResolutionReference = (25, 25, 25);    

    #Orientation of the Data set wrt reference 
    #(-axis will invert the orientation, for other hemisphere use (-1, 2, 3), to exchnge x,y use (2,1,3) etc)
    Orientation = None;
    


##############################################################################
# Analysis
##############################################################################

class VoxelizationParameter(object):
    """Parameter to calculate density voxelization"""
    
    # Define bounds of the volume to be voxelized
    MinX = 0;
    MaxX = 0;
    MinY = 0;
    MaxY = 0;    
    MinZ = 0;
    MaxZ = 0;
    Distance = 200;  #Radius of the sphere
    
    #Spacing
    XSpacing = 0;
    YSpacing = 0;
    ZSpacing = 0;
    
    #Image dimensions
    NX = 0;
    NY = 0;
    NZ = 0;
    
    Type = 'Spherical' ; # Spherical, , 'Rectangular, Gaussian'






##############################################################################
# Full Parameter
##############################################################################

class Parameter(object):
    """ Parameter for full processing chain"""
    
    #Input files
    DataSource = DataSourceParameter();
    
    #Elastix alignment parameter
    Alignment = AlignmentParameter();
    
    #Resampling 
    Resampling = ResamplingParameter();
    
    #Image Rpcoessing parameter
    ImageProcessing = ImageProcessingParameter();
    
    #ParallelProcessing
    StackProcessing = StackProcessingParameter();
    
    #Voxelization
    Voxelization = VoxelizationParameter();
    
    #Plot intermediate results for debugging
    Verbose = False;
    


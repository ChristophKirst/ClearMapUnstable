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
    fn = os.path.split(self.__file__);
    return fn[0];



class DataSourceParameter(object):
    """Parameter for files"""
    
    #File name or file pattern of raw image data
    ImageFile = None;
         
    #Co-ordinate ranges of raw data
    ZRange = all;
    XRange = all;
    YRange = all;
        

class SpotDetectionParameter(object):
    
    # Background correctoin: None or (x,y) which is size of disk for gray scale opening
    Background = (15,15);
    
    # Spot Detection via Difference of Gaussians (DoG) filter: (x,y,z) size
    Dog = (11, 7, 7);
    
    #h of h-max transform
    HMax = 20;
    
    #Threshold for min intensity at center to be counted as cell (should be similar to the h max)
    Threshold = 20;
    
    
class IlastikParameter(object):
    
    #ilastik 0.5 path
    Ilastik = None;
    
    #ilastic classifier to use
    Classifier = os.path.join(self.iDISCOPath(), '/Test/Ilastik/classifier.h5');  
    
    # Rescaling of images
    Rescale = None;
    
    # Background correctoin: None or (x,y) which is size of disk for gray scale opening
    Background = (15,15);




class ImageProcessingParameter(object):
    """ Parameter for image processing chain"""
    
    Method = 'SpotDetection';
    
    Parameter = SpotDetectionParameter();
    
    #File name for cell coordinates csv, vtk or  ims extension
    CellCoordinateFile = None;



class AlignmentParameter(object):
    
    #directory of the alignment result
    AlignmentDirectory = None;
    
    #Elastix binary
    Elastix = '/usr/local/elastix/bin/elastix'
    
    #Transformix binary
    Transformix = '/usr/local/elastix/bin/transformix'
    
    Lib = '/usr/local/elastix/lib'
    
    #moving and reference images
    MovingImage = os.path.join(self.iDISCOPath(), '/Test/Data/Elastix/150524_0_8X-s3-20HFautofluor_18-51-1-warpable.tif');
    FixedImage  = os.path.join(self.iDISCOPath(), '/Test/Data/Elastix/OstenRefARA_v2_lowerHalf.tif');
    FixedImageMask = None;
    
    #elastix parameter files for alignment
    AffineParameterFile  = os.path.join(self.iDISCOPath(), '/Test/Elastix/ElastixParamterAffine.txt');
    BSplineParameterFile = os.path.join(self.iDISCOPath(), '/Test/Elastix/ElastixParamterBSpline.txt');
    
    
    
class StackProcessingParameter(object):
    
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
    



class Parameter(object):
    """ Parameter for full processing chain"""
    
    #Image Rpcoessing parameter
    ImageProcessing = ImageProcessingParameter();
    
    #ParallelProcessing
    StackProcessing = StackProcessingParameter();
    
    #Input files
    DataSource = DataSourceParameter();
    
    #Elastix alignment parameter
    Alignment = AlignmentParameter();
    
    #Plot intermediate results for debugging
    Verbose = False;
    


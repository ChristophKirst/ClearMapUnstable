# -*- coding: utf-8 -*-

"""
Parameter 

Notes:
    - defines standard/default parameters

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import os
import iDISCO.Settings as settings

#from iDISCO.Utils.ParameterTools import joinParameter

##############################################################################
# Image Processing
##############################################################################

## Paramters for cell detection using spot detection alorithm 
SpotDetectionParameter = {
    # Background correctoin: None or (y,x) which is size of disk for gray scale opening
    "backgroundSize" : (15,15),
    
    # Spot Detection via Difference of Gaussians (DoG) filter: (y,x,z) size
    "dogSize" : (7, 7, 11),
    
    #h of h-max transform
    "hMax" : 20,
    
    #intensity detection   
    "intensityMethod"  : 'Max',  #None -> intensity of pixel of center, alternatively string of numpy array method that returns a single number
    "intensitySize"    : (3,3,3),  # size of box in (y,x,z) to include in intensity determination
    
    #Threshold for min intensity at center to be counted as cell (should be similar to the h max)
    "threshold" : 20
    };
    

## Paramters for cell detection using Ilastik classification 
IlastikParameter = {
    #ilastic classifier to use
    "classifier" : os.path.join(settings.IDISCOPath, '/Test/Ilastik/classifier.h5'),
    
    # Rescaling of images
    "rescale" : None,
    
    # Background correctoin: None or (y,x) which is size of disk for gray scale opening
    "backgroundSize" : (15,15)
    };


## Parameter for processing a stack in parallel
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 2,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "cChunkOverlap" : 15,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : True,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
    "chunkOptimizationSize" : all
    };


##############################################################################
# Alignment / Morphin / Resmapling
##############################################################################

## Parameter for Elastix alignment
AlignmentParameter = {
    #directory of the alignment result
    "alignmentDirectory" : None,
            
    #moving and reference images
    "movingImage" : os.path.join(settings.IDISCOPath, '/Test/Data/Elastix/150524_0_8X-s3-20HFautofluor_18-51-1-warpable.tif'),
    "fixedImage"  : os.path.join(settings.IDISCOPath, '/Test/Data/Elastix/OstenRefARA_v2_lowerHalf.tif'),
    "fixedImageMask" : None,
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(settings.IDISCOPath, '/Test/Elastix/ElastixParameterAffine.txt'),
    "bSplineParameterFile" : os.path.join(settings.IDISCOPath, '/Test/Elastix/ElastixParameterBSpline.txt'),
    };


## Parameter for resampling data
ResamplingParameter = {
    
    #Data source and output file
    "source" : None,
    "sink"   : None,
    
    #Resolution of the raw data (in um / pixel) as (x,y,z)
    "resolutionSource" : (4.0625, 4.0625, 3),

    #Resolution of the Reference / Atlas (in um/ pixel) as (x,y,z)
    "resolutionSink" : (25, 25, 25),

    #Orientation of the Data set wrt reference as (x=1,y=2,z=3)
    #(-axis will invert the orientation, for other hemisphere use (-1, 2, 3), to exchnge x,y use (2,1,3) etc)
    "orientation" : None
    };



##############################################################################
# Analysis
##############################################################################

## Parameter to calculate density voxelization
VoxelizationParameter = {
    #Method to voxelize
    "method" : 'Spherical', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized
    "voxelizationSize" : (1,1,1),  
    };




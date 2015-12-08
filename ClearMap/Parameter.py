# -*- coding: utf-8 -*-
"""
ClearMap default parameter module.

This module defines default parameter used by various sub-packages.

See Also:
   :mod:`~ClearMap.Settings`

Author
""""""
   Christoph Kirst, The Rockefeller University, 2015
"""

import os
import ClearMap.Settings as settings

#from ClearMap.Utils.ParameterTools import joinParameter

##############################################################################
# Image Processing
##############################################################################

## Paramters for cell detection using spot detection algorithm 
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
"""
dict: Paramters for cell detection using the spot detection algorithm

   * "backgroundSize" : Background correctoin: None or (y,x) which is size of disk for gray scale opening
    
   * "dogSize" : Spot Detection via Difference of Gaussians (DoG) filter: (y,x,z) size

   * "hMax" : h of h-max transform

   * "intensityMethod"  : None -> intensity of pixel of center, alternatively string of numpy array method that returns a single number
   * "intensitySize"    : size of box in (y,x,z) to include in intensity determination    
  
   * "threshold" : Threshold for min intensity at center to be counted as cell (should be similar to the h max)
   
See Also:
   :const:`IlastikParameter`, :const:`StackProcessingParameter`
"""     

## Paramters for cell detection using Ilastik classification 
IlastikParameter = {
    #ilastic classifier to use
    "classifier" : os.path.join(settings.ClearMapPath, '/Test/Ilastik/classifier.h5'),
    
    # Rescaling of images
    "rescale" : None,
    
    # Background correctoin: None or (y,x) which is size of disk for gray scale opening
    "backgroundSize" : (15,15)
    };
"""
dict: Paramters for cell detection using Ilastik classification

   * "classifier": ilastic classifier to use
     
   * "rescale": rescale images before classification
    
   * "backgroundSize": Background correctoin: None or (y,x) which is size of disk for gray scale opening

See Also:
   :const:`SpotDetectionParameter`, :const:`StackProcessingParameter`
"""

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
"""
dict: Parameter for processing an image stack in parallel
    *  "processes": max number of parallel processes

    * "chunkSizeMax" : maximal chunk size in z 
    * "chunkSizeMin" : minimal chunk size in z,
    * "chunkOverlap" : overlap between two chunks,

    * "chunkOptimization": optimize chunk size and number to number of processes
    * "chunkOptimizationSize": increase chunk size for optimizaition (True, False or all = automatic)

See Also:
   :const:`SpotDetectionParameter`, :const:`IlastikParameter`
"""

##############################################################################
# Alignment / Morphin / Resmapling
##############################################################################

## Parameter for Elastix alignment
AlignmentParameter = {
    #directory of the alignment result
    "alignmentDirectory" : None,
            
    #moving and reference images
    "movingImage" : os.path.join(settings.ClearMapPath, '/Test/Data/Elastix/150524_0_8X-s3-20HFautofluor_18-51-1-warpable.tif'),
    "fixedImage"  : os.path.join(settings.ClearMapPath, '/Test/Data/Elastix/OstenRefARA_v2_lowerHalf.tif'),
    "fixedImageMask" : None,
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(settings.ClearMapPath, '/Test/Elastix/ElastixParameterAffine.txt'),
    "bSplineParameterFile" : os.path.join(settings.ClearMapPath, '/Test/Elastix/ElastixParameterBSpline.txt'),
    };
"""
dict: Parameter for Elastix alignment

  * "alignmentDirectory" : directory to save the alignment result
              
  * "movingImage": image to be aligned 
  * "fixedImage":  reference image
  
  * "affineParameterFile": elastix parameter files for affine alignment
  * "bSplineParameterFile" : elastix parameter files for non-linear alignment

See Also:
   :mod:`~ClearMap.Alignment.Elastix`
"""

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
"""
dict: Parameter for resampling data

    * "source" : data source file
    * "sink"   : data output file
    
    "resolutionSource": resolution of the raw data (in um / pixel) as (x,y,z)

    "resolutionSink" : resolution of the reference / atlas image (in um/ pixel) as (x,y,z)

    "orientation" : Orientation of the data set wrt reference as (x=1,y=2,z=3) 
     (-axis will invert the orientation, for other hemisphere use (-1, 2, 3), to exchnge x,y use (2,1,3) etc)
  
See Also:
   :mod:`~ClearMap.Alignment.Resampling`
"""


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
"""
dict: Parameter to calculate density voxelization

    * "method": Method to voxelize: 'Spherical','Rectangular, 'Gaussian'
    * "voxelizationSize": max size of the volume to be voxelized
    
See Also:
   :mod:`~ClearMap.Analysis.voxelization`
"""




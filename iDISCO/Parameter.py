# -*- coding: utf-8 -*-

"""
Example Parameter file

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

class IOParameter(object):
    
    #File name or file pattern of raw image data
    ImageFile = None;
    
    #directory of the alignment result
    AlignmentDirectory = None;

    #File name for cell coordinates csv, vtk or  ims extension
    CellCoordinateFile = None;
    

class ImageProcessingParameter(object):
    """ Parameter for image processing chain"""
    
    # Background correctoin: None or (x,y) which is size of disk for gray scale opening
    Background = (15,15);
    
    # Spot Detection via Difference of Gaussians (DoG) filter: (x,y,z) size
    Dog = (11, 7, 7);
    
    #h of h-max transform
    HMax = 20;
    
    #Threshold for min intensity at center to be counted as cell (should be similar to the h max)
    Threshold = 20;
    
    #Ilastik classifier
    IlastikClassifier = None;
    


class AlignmentParameter(object):
    
    #Elastix binary
    Elastix = '/usr/local/elastix/bin/elastix'
    
    #TRansformix binary
    Transformix = '/usr/local/elastix/bin/transformix'
    
    #moving and reference images
    MovingImage = './Test/Data/Elastix/150524_0_8X-s3-20HFautofluor_18-51-1-warpable.tif'
    FixedImage  = './Test/Data/Elastix/OstenRefARA_v2_lowerHalf.tif'
    FixedImageMask = None;
    
    #elastix parameter files for alignment
    AffineParameterFile  = './Test/Elastix/ElastixParamterAffine.txt'
    BSplineParameterFile = './Test/Elastix/ElastixParamterBSpline.txt'
    

  
class Parameter(object):
    """ Parameter for full processing chain"""
    
    #Image Rpcoessing parameter
    ImageProcessing = ImageProcessingParameter();
    
    #Input files
    IO = IOParameter();
    
    #Elastix alignment parameter
    Alignment = AlignmentParameter();
    
    #Plot intermediate results for debugging
    Verbose = False;
    


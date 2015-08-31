# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 17:25:40 2015

@author: mtllab
"""

import iDISCO.Alignment.Elastix as elx
import iDISCO.Alignment.Resampling as rsp
import os, numpy

import iDISCO.IO.IO as io

BaseDirectory = '/home/mtllab/Documents/MRI';

#Path to registration parameters and atlases
PathReg        = '/home/mtllab/Documents/warping';



AffineAlignmentParameter = {            
    #moving and reference images
    "movingImage" : os.path.join(BaseDirectory, 'lightsheet.tif'),
    "fixedImage"  :  os.path.join(BaseDirectory, 'MRI.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(PathReg, 'Par0000affine.txt'),
    "bSplineParameterFile" : None,
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_affine')
    }; 

resultDirectory  = elx.alignData(**AffineAlignmentParameter);

lightsheetTransformed = elx.transformData(AffineAlignmentParameter["movingImage"], sink =  os.path.join(BaseDirectory, 'lightsheet_affine.tif'), transformDirectory = resultDirectory, resultDirectory = None)



NonLinearAlignmentParameter = {            
    #moving and reference images
    "movingImage" :  lightsheetTransformed,
    "fixedImage"  :  os.path.join(BaseDirectory, 'MRI.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : None,
    "bSplineParameterFile" : os.path.join(PathReg, 'Par0000bspline.txt'),
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_bspline')
    }; 

resultDirectory  = elx.alignData(**NonLinearAlignmentParameter);




elx.deformationField(transformDirectory = NonLinearAlignmentParameter["resultDirectory"], resultDirectory = '/home/mtllab/Documents/MRI/deformationField')

deformation = io.readData('/home/mtllab/Documents/MRI/deformationField/deformationField.mhd');

distance = elx.deformationDistance(deformation);

io.writeData('/home/mtllab/Documents/MRI/deformation.tif', deformation)

io.writeData('/home/mtllab/Documents/MRI/distance.mhd', distance)



tf = elx.getTransformParameterFile(AffineAlignmentParameter["resultDirectory"]);


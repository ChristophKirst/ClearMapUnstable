# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 17:25:40 2015

@author: mtllab
"""

import iDISCO.Alignment.Elastix as elx
import iDISCO.Alignment.Resampling as rsp
import os, numpy
from iDISCO.Alignment.Resampling import resampleData;
import iDISCO.IO.IO as io

BaseDirectory = '/home/mtllab/Documents/MRI';

#Path to registration parameters and atlases
PathReg        = '/home/mtllab/Documents/warping';


ResamplingParameter = { "processes": 12,
                       "source" : '/home/mtllab/Documents/MRI/lightsheet/23/151025_0_8X_autofluo_16-04-39/16-04-39_0_8X_autofluo_UltraII_C00_xyz-Table Z\d{4}.ome.tif',
                       "resolutionSource" : (4.0625, 4.0625, 3),
                       "resolutionSink"   : (25, 25, 25),
                       "orientation" : (3, 1, 2),
                       "sink"   : os.path.join(BaseDirectory, 'LSFM_23_test.tif')};

resampleData(**ResamplingParameter);

##########################################################3

AffineAlignmentParameter = {            
    #moving and reference images
    "movingImage" : os.path.join(BaseDirectory, 'lightsheet.tif'),
    "fixedImage"  :  os.path.join(BaseDirectory, 'MRI.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(PathReg, 'Par0000affine.txt'),
    "bSplineParameterFile" : None,
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_affine_WB')
    }; 

resultDirectory  = elx.alignData(**AffineAlignmentParameter);

lightsheetTransformed = elx.transformData(AffineAlignmentParameter["movingImage"], sink =  os.path.join(BaseDirectory, 'LSFM_WB_affine.tif'), transformDirectory = resultDirectory, resultDirectory = None)

NonLinearAlignmentParameter = {            
    #moving and reference images
    "movingImage" :  lightsheetTransformed,
    "fixedImage"  :  os.path.join(BaseDirectory, 'mri19.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : None,
    "bSplineParameterFile" : os.path.join(PathReg, 'Par0000bspline.txt'),
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_bspline_19')
    }; 

resultDirectory  = elx.alignData(**NonLinearAlignmentParameter);

elx.deformationField(transformDirectory = NonLinearAlignmentParameter["resultDirectory"], resultDirectory = os.path.join(BaseDirectory, 'deformationfield_19'))

deformation = io.readData(os.path.join(BaseDirectory, 'deformationfield_19/deformationField.mhd'));

distance = elx.deformationDistance(deformation);

io.writeData('/home/mtllab/Documents/MRI/deformation_19.tif', deformation)

io.writeData('/home/mtllab/Documents/MRI/distance_19.mhd', distance)
########################################################################################

##########################################################3

AffineAlignmentParameter = {            
    #moving and reference images
    "movingImage" : os.path.join(BaseDirectory, 'LSFM_20_resliced.tif'),
    "fixedImage"  :  os.path.join(BaseDirectory, 'mri20.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(PathReg, 'Par0000affine.txt'),
    "bSplineParameterFile" : None,
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_affine_20')
    }; 

resultDirectory  = elx.alignData(**AffineAlignmentParameter);

lightsheetTransformed = elx.transformData(AffineAlignmentParameter["movingImage"], sink =  os.path.join(BaseDirectory, 'LSFM_20_affine.tif'), transformDirectory = resultDirectory, resultDirectory = None)

NonLinearAlignmentParameter = {            
    #moving and reference images
    "movingImage" :  lightsheetTransformed,
    "fixedImage"  :  os.path.join(BaseDirectory, 'mri20.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : None,
    "bSplineParameterFile" : os.path.join(PathReg, 'Par0000bspline.txt'),
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_bspline_20')
    }; 

resultDirectory  = elx.alignData(**NonLinearAlignmentParameter);

elx.deformationField(transformDirectory = NonLinearAlignmentParameter["resultDirectory"], resultDirectory = os.path.join(BaseDirectory, 'deformationfield_20'))

deformation = io.readData(os.path.join(BaseDirectory, 'deformationfield_20/deformationField.mhd'));

distance = elx.deformationDistance(deformation);

io.writeData('/home/mtllab/Documents/MRI/deformation_20.tif', deformation)

io.writeData('/home/mtllab/Documents/MRI/distance_20.mhd', distance)
########################################################################################
##########################################################3

AffineAlignmentParameter = {            
    #moving and reference images
    "movingImage" : os.path.join(BaseDirectory, 'LSFM_21_resliced.tif'),
    "fixedImage"  :  os.path.join(BaseDirectory, 'mri21.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(PathReg, 'Par0000affine.txt'),
    "bSplineParameterFile" : None,
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_affine_21')
    }; 

resultDirectory  = elx.alignData(**AffineAlignmentParameter);

lightsheetTransformed = elx.transformData(AffineAlignmentParameter["movingImage"], sink =  os.path.join(BaseDirectory, 'LSFM_21_affine.tif'), transformDirectory = resultDirectory, resultDirectory = None)

NonLinearAlignmentParameter = {            
    #moving and reference images
    "movingImage" :  lightsheetTransformed,
    "fixedImage"  :  os.path.join(BaseDirectory, 'mri21.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : None,
    "bSplineParameterFile" : os.path.join(PathReg, 'Par0000bspline.txt'),
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_bspline_21')
    }; 

resultDirectory  = elx.alignData(**NonLinearAlignmentParameter);

elx.deformationField(transformDirectory = NonLinearAlignmentParameter["resultDirectory"], resultDirectory = os.path.join(BaseDirectory, 'deformationfield_21'))

deformation = io.readData(os.path.join(BaseDirectory, 'deformationfield_21/deformationField.mhd'));

distance = elx.deformationDistance(deformation);

io.writeData('/home/mtllab/Documents/MRI/deformation_21.tif', deformation)

io.writeData('/home/mtllab/Documents/MRI/distance_21.mhd', distance)
########################################################################################
##########################################################3

AffineAlignmentParameter = {            
    #moving and reference images
    "movingImage" : os.path.join(BaseDirectory, 'LSFM_22_resliced.tif'),
    "fixedImage"  :  os.path.join(BaseDirectory, 'mri22.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(PathReg, 'Par0000affine.txt'),
    "bSplineParameterFile" : None,
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_affine_22')
    }; 

resultDirectory  = elx.alignData(**AffineAlignmentParameter);

lightsheetTransformed = elx.transformData(AffineAlignmentParameter["movingImage"], sink =  os.path.join(BaseDirectory, 'LSFM_22_affine.tif'), transformDirectory = resultDirectory, resultDirectory = None)

NonLinearAlignmentParameter = {            
    #moving and reference images
    "movingImage" :  lightsheetTransformed,
    "fixedImage"  :  os.path.join(BaseDirectory, 'mri22.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : None,
    "bSplineParameterFile" : os.path.join(PathReg, 'Par0000bspline.txt'),
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_bspline_22')
    }; 

resultDirectory  = elx.alignData(**NonLinearAlignmentParameter);

elx.deformationField(transformDirectory = NonLinearAlignmentParameter["resultDirectory"], resultDirectory = os.path.join(BaseDirectory, 'deformationfield_22'))

deformation = io.readData(os.path.join(BaseDirectory, 'deformationfield_22/deformationField.mhd'));

distance = elx.deformationDistance(deformation);

io.writeData('/home/mtllab/Documents/MRI/deformation_22.tif', deformation)

io.writeData('/home/mtllab/Documents/MRI/distance_22.mhd', distance)
########################################################################################
##########################################################3

AffineAlignmentParameter = {            
    #moving and reference images
    "movingImage" : os.path.join(BaseDirectory, 'LSFM_23_resliced.tif'),
    "fixedImage"  :  os.path.join(BaseDirectory, 'mri23.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : os.path.join(PathReg, 'Par0000affine.txt'),
    "bSplineParameterFile" : None,
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_affine_23')
    }; 

resultDirectory  = elx.alignData(**AffineAlignmentParameter);

lightsheetTransformed = elx.transformData(AffineAlignmentParameter["movingImage"], sink =  os.path.join(BaseDirectory, 'LSFM_23_affine.tif'), transformDirectory = resultDirectory, resultDirectory = None)

NonLinearAlignmentParameter = {            
    #moving and reference images
    "movingImage" :  lightsheetTransformed,
    "fixedImage"  :  os.path.join(BaseDirectory, 'mri23.tif'),
    
    #elastix parameter files for alignment
    "affineParameterFile"  : None,
    "bSplineParameterFile" : os.path.join(PathReg, 'Par0000bspline.txt'),
    
    #directory of the alignment result
    "resultDirectory" :  os.path.join(BaseDirectory, 'elastix_bspline_23')
    }; 

resultDirectory  = elx.alignData(**NonLinearAlignmentParameter);

elx.deformationField(transformDirectory = NonLinearAlignmentParameter["resultDirectory"], resultDirectory = os.path.join(BaseDirectory, 'deformationfield_23'))

deformation = io.readData(os.path.join(BaseDirectory, 'deformationfield_23/deformationField.mhd'));

distance = elx.deformationDistance(deformation);

io.writeData('/home/mtllab/Documents/MRI/deformation_23.tif', deformation)

io.writeData('/home/mtllab/Documents/MRI/distance_23.mhd', distance)
########################################################################################

#19: x=[[0.902552, 0.076987, -0.170393], [-0.139207, 0.994968, -0.069795], [0.213558, 0.091294, 0.925352]]
#20:  x=[[0.908789, 0.017137, -0.141930],[-0.059510, 0.982788, -0.101138],[0.139824, 0.124675, 0.943861]]
#21: z=[[0.919094, -0.012410, -0.295410],[ 0.040059, 0.992951, 0.037802],[ 0.278923 ,-0.060937, 0.994752]]
#22: v=[[0.889665, 0.197466, -0.180494],[ -0.295030, 1.035758, -0.015283],[ 0.187451, 0.042624, 0.858066]]
#23: x=[[0.881118, 0.182486, -0.258064],[ -0.334966 ,1.060279, -0.092897 ],[0.242058, 0.159615, 0.725423]]


#tf = elx.getTransformParameterFile(AffineAlignmentParameter["resultDirectory"]);

#################################################3

import matplotlib.pyplot as matpp
distance = io.readData('/home/mtllab/Documents/MRI/distance_23.mhd');
d = distance[60:200,60:210,150:270];
matpp.hist(d.flatten(), bins=1000, histtype="step");
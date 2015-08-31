# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 15:40:34 2015

@author: mtllab
"""

import iDISCO.Alignment.Elastix as elx
import iDISCO.Alignment.Resampling as rsp
import iDISCO.Analysis.Label as lbl
import os, numpy



#Path to registration parameters and atlases
PathReg        = '/home/mtllab/Documents/warping';
AtlasFile      = os.path.join(PathReg, 'half_template_25_right.tif');
AnnotationFile = os.path.join(PathReg, 'annotation_25_right.tif');

baseDirectory= '/home/mtllab/Documents/Haloperidol/1strun/Haloperidol'



resultFile = elx.transformData(AnnotationFile, sink = [], 
                  transformDirectory = os.path.join(baseDirectory, 'elastix_auto_to_atlas'), resultDirectory = os.path.join(baseDirectory, 'elastix_atlas_to_auto'))


resultFile = os.path.join(baseDirectory, 'elastix_atlas_to_auto/result.mhd');

resultFile2 = rsp.resampleData(resultFile, sink = os.path.join(baseDirectory, 'atlas_to_cfos_resampled.tif'),  orientation = None, 
                 dataSizeSink = os.path.join(baseDirectory, 'cfos_resampled.tif'), interpolation = None, processes = 12)


resultFile2 = os.path.join(baseDirectory, 'atlas_to_cfos_resampled.tif')

resultFile3 = elx.transformData(resultFile2, sink = [], 
                  transformDirectory = os.path.join(baseDirectory, 'elastix_cfos_to_auto'), resultDirectory = os.path.join(baseDirectory, 'elastix_atlas_to_cfos_resampled'));



resultFile3 = os.path.join(baseDirectory, 'elastix_atlas_to_cfos_resampled/result.mhd');


resultFile4 = rsp.resampleDataInverse(resultFile3, source = os.path.join(baseDirectory, 'atlas_to_cfos/Z\d{4}.tif'),  orientation = None, 
                 dataSizeSource = os.path.join(baseDirectory, '150630_0_8X-cfos_18-32-28/18-32-28_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif'), interpolation = None, processes = 12)


import iDISCO.IO.IO as io

io.dataSize(os.path.join(baseDirectory, '150630_0_8X-cfos_18-32-28/18-32-28_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif'))
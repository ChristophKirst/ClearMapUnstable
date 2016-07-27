# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 13:48:09 2016

@author: ckirst
"""

# test cropping


import ClearMap.IO as io
reload(io);


import ClearMap.Settings as settings

import ClearMap.Alignment.Stitching as st

import os

import numpy as np


datadir ='/home/mtllab/Documents/th/';

fn = os.path.join(datadir, r'160412_mosaic_15-20-19/15-20-19_mosaic_UltraII\[(?P<row>\d{2}) x (?P<col>\d{2})\]_C00_xyz-Table Z(?P<z>\d{4}).ome.tif')

_, gr = st.findFileList(fn , sort = True, groups = ['row','col'], absolute = True)
groups = [];
for i in range(gr.shape[1]):
    groups.append(np.unique(gr[:,i]));

print groups

for i in groups[0]:
    for j in groups[1]:
        fileExpression = os.path.join(datadir, r'160412_mosaic_15-20-19/15-20-19_mosaic_UltraII\[%s x %s]_C00_xyz-Table Z\d{4}.ome.tif' % (i,j))

        io.dataSize(fileExpression)

        io.readMetaData(fileExpression, info = ['size', 'overlap', 'resolution'])


        import ClearMap.IO.FileList as fl;
        reload(fl)

        fncrop = os.path.join(datadir, r'cropped/15-20-19_mosaic_UltraII_%s_x_%s_C00_xyz-Table Z\d{4}.ome.tif' % (i,j))


        fc = fl.cropData(fileExpression, fncrop, x = (400, -400), y = (550, -550),  adjustOverlap = True, processes = all)
        #fc1 = fl.firstFile(fc)
        #io.readMetaData(fc1, info = ['overlap', 'resolution', 'size']);




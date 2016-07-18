# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 13:48:09 2016

@author: ckirst
"""

# test cropping


import ClearMap.IO as io
reload(io);


import ClearMap.Settings as settings

import os

fn = os.path.join(settings.ClearMapPath, 'Test/Data/Stitching/14-22-21_TH-s3-ov45-na004_UltraII\[00 x 00\]_C00_xyz-Table Z\d{4}.ome.tif');

io.dataSize(fn)

io.readMetaData(fn, info = ['size', 'overlap', 'resolution'])


import ClearMap.IO.FileList as fl;
reload(fl)

fncrop = os.path.join(settings.ClearMapPath, 'Test/Data/Cropped/crop_test_\d{4}.tif');

fc = fl.cropData(fn, fncrop, x = (500, -500), y = (500, -500),  adjustOverlap = True, processes = all)
fc1 = fl.firstFile(fc)
io.readMetaData(fc1, info = ['overlap', 'resolution', 'size']);

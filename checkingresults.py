# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 16:03:57 2015

@author: mtllab
"""

import iDISCO.Visualization.Plot as plt;
import iDISCO.IO.IO as io;

dataSource =  '/home/mtllab/Documents/Haloperidol/1266/150818_0_8X-fos_14-44-03/14-44-03_0_8X-fos_UltraII_C00_xyz-Table Z\d{4}.ome.tif'
pointSource= '/home/mtllab/Documents/Haloperidol/1266/cells-allpoints.npy';
x = (715,1100);
y = (978,1371);
z = (667,742);


data = plt.overlayPoints(dataSource, pointSource, x = x, y = y, z = z, pointColor = None);
io.writeData('/home/mtllab/Documents/Haloperidol/1266/checkpoints.tif', data)

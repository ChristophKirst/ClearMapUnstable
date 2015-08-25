# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 16:03:57 2015

@author: mtllab
"""

import iDISCO.Visualization.Plot as plt;
import iDISCO.IO.IO as io;

dataSource =  '/home/mtllab/Documents/Haloperidol/1267/150818_0_8X-cfos_17-14-34/17-14-34_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif'
pointSource= '/home/mtllab/Documents/Haloperidol/1267test/cells.csv';
x = (815,1000);
y = (1078,1271);
z = (667,742);



p, i = io.readPoints(('/home/mtllab/Documents/Haloperidol/1267test/cells.csv', '/home/mtllab/Documents/Haloperidol/1267test/intensities.csv'))
p2, i2 = io.pointsToRange((p,i), dataSource, x = x, y = y, z = z)
p3 = p2 -io.pointShiftFromRange(dataSource, x = x, y = y, z = z)
p3a = p3.astype(int);


p3
p3a
i2


p4, i4 = io.pointsToRange((p3,i2), dataSource, x = (5, 80))
p4
i4



data = plt.overlayPoints(dataSource, pointSource, x = x, y = y, z = z, pointColor = None);
io.writeData('/home/mtllab/Documents/Haloperidol/1267/check.tif', data)

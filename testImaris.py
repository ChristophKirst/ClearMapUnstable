# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 20:11:49 2015

@author: nicolas
"""

import iDISCO.IO.IO as io

import iDISCO.IO.Imaris as ims

import numpy


points = io.readPoints('/home/nicolas/Windows/Nico/cfosRegistrations/haloperidol/cells.npy');

imsFile = '/home/nicolas/Windows/Nico/cfosRegistrations/haloperidol/1266_small.ome.ims';


print io.dataSize(imsFile)

data = io.readData(imsFile, resolution = 0);
print data.shape

points = numpy.array([[0,0,0], [234, 497, 472], [10,20,30], [40,90, 120]]);

#points = numpy.array([[0,0,0], [951, 2019, 1416], [100,200,300]]);

io.writePoints(imsFile, points)






h5file = ims.openFile(imsFile)

h5file.close()




points = io.readPoints('/home/nicolas/Windows/Nico/cfosRegistrations/haloperidol/cells.npy');

imsFile = '/home/nicolas/Windows/Nico/cfosRegistrations/haloperidol/1266.ome.ims';

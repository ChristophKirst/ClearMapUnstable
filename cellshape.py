# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 22:41:09 2015

@author: mtllab
"""

points = io.readPoints('/home/mtllab/Documents/Haloperidol/1266/cells-allpoints.npy');

r = {'x' : (1000,1500), 'y' : (700,1200), 'z' : (500,800)};

sh = io.pointShiftFromRange(cFosFile, **cFosFileRange);

points = points - sh;

import iDISCO.Analysis.Voxelization as vox


img = vox.voxelizePixel(points, io.dataSize(cFosFile, **r))

io.writeData('/home/mtllab/Documents/Haloperidol/1266/cells.tif', img)



ic = io.readPoints('/home/mtllab/Documents/Haloperidol/1266/intensities-allpoints.npy')

import matplotlib.pylab as mpl


mpl.hist(ic[:,3], bins = 100, range = (0, 100))


iid = ic[:,3] > 30;

iid2 = ic[:,3] > 2000;
numpy.sum(iid2)

ii = ic[iid,0];


mpl.figure()
mpl.hist(ii, bins = 100)


numpy.sum(iid)
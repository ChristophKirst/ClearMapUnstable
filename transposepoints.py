# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 15:42:29 2015

@author: mtllab
"""

import iDISCO.IO.IO as io;

points = io.readPoints('/home/mtllab/Documents/Haloperidol/1267/cells.csv');
newpoints = points[:,(1,0,2)];
io.writePoints('/home/mtllab/Documents/Haloperidol/1267/cells-transposed.csv', newpoints);


points = io.readPoints('/home/mtllab/Documents/Haloperidol/1268/cells.csv');
newpoints = points[:,(1,0,2)];
io.writePoints('/home/mtllab/Documents/Haloperidol/1268/cells-transposed.csv', newpoints);


points = io.readPoints('/home/mtllab/Documents/Haloperidol/1269/cells.csv');
newpoints = points[:,(1,0,2)];
io.writePoints('/home/mtllab/Documents/Haloperidol/1269/cells-transposed.csv', newpoints);


points = io.readPoints('/home/mtllab/Documents/Haloperidol/1270/cells.csv');
newpoints = points[:,(1,0,2)];
io.writePoints('/home/mtllab/Documents/Haloperidol/1270/cells-transposed.csv', newpoints);


points = io.readPoints('/home/mtllab/Documents/Haloperidol/1271/cells.csv');
newpoints = points[:,(1,0,2)];
io.writePoints('/home/mtllab/Documents/Haloperidol/1271/cells-transposed.csv', newpoints);


points = io.readPoints('/home/mtllab/Documents/Haloperidol/1272/cells.csv');
newpoints = points[:,(1,0,2)];
io.writePoints('/home/mtllab/Documents/Haloperidol/1272/cells-transposed.csv', newpoints);


points = io.readPoints('/home/mtllab/Documents/Haloperidol/1273/cells.csv');
newpoints = points[:,(1,0,2)];
io.writePoints('/home/mtllab/Documents/Haloperidol/1273/cells-transposed.csv', newpoints);


points = io.readPoints('/home/mtllab/Documents/Haloperidol/1274/cells.csv');
newpoints = points[:,(1,0,2)];
io.writePoints('/home/mtllab/Documents/Haloperidol/1274/cells-transposed.csv', newpoints);


points = io.readPoints('/home/mtllab/Documents/Haloperidol/1275/cells.csv');
newpoints = points[:,(1,0,2)];
io.writePoints('/home/mtllab/Documents/Haloperidol/1275/cells-transposed.csv', newpoints);


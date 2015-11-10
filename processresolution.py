# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 20:23:41 2015

@author: mtllab
"""


# -*- coding: utf-8 -*-
"""
Created on Sat Aug  8 14:58:07 2015

@author: mtllab"""

import iDISCO.Visualization.Plot as plt;


execfile('0_8Xresolution.py')

detectCells(**ImageProcessingParameter);

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints-0_8X.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints-0_8X.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (60, 500), row = (3,3));
io.writePoints((os.path.join(BaseDirectory, 'cells-0_8X.npy'), os.path.join(BaseDirectory,  'intensities-0_8X.npy')), (points, intensities));

#################### Heat map generation

data = plt.overlayPoints('/home/mtllab/Documents/countingaccuracy/resolution/151012_0_8Xs3_16-09-57.tif', points, x = all, y = all, z = all, pointColor = None);
io.writeData(os.path.join(BaseDirectory, 'cells0_8X_points0_8X.tif'), data)


#
#
######################
######################
######################
######################
######################
######################
#

execfile('2Xresolution.py')

detectCells(**ImageProcessingParameter);

points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints-2X.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints-2X.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (350, 2850), row = (3,3));
io.writePoints((os.path.join(BaseDirectory, 'cells-2X.npy'), os.path.join(BaseDirectory,  'intensities-2X.npy')), (points, intensities));

#################### Heat map generation

ResamplingPointsParameter = { "processes": 12 };
ResamplingPointsParameter["pointSource"] = os.path.join(BaseDirectory, 'cells-2X.npy');

ResamplingPointsParameter["resolutionSource"] = (1.63,1.63,3);
ResamplingPointsParameter["resolutionSink"]   = (4.065,4.065,3);

ResamplingPointsParameter["dataSizeSource"] = '/home/mtllab/Documents/countingaccuracy/resolution/151012_2Xs3_17-06-22.tif';
ResamplingPointsParameter["pointSink"]  = None;



points = resamplePoints(**ResamplingPointsParameter);


data = plt.overlayPoints('/home/mtllab/Documents/countingaccuracy/resolution/151012_2Xs3_17-06-22.tif', points, x = all, y = all, z = all, pointColor = None);
io.writeData(os.path.join(BaseDirectory, 'cells2X_points2X.tif'), data)
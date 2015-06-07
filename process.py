# -*- coding: utf-8 -*-
"""
Processing Brain Samples - Main Example Script


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import numpy as np
import matplotlib as mpl

import iDISCO.Visualization.Plot as iplt
#import iDISCO.Segmentation.SpotDetection as iseg
import iDISCO.IO.Imaris as io


# open data
fn = '/home/ckirst/Science/Projects/BrainActivityMap/Data/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims';
fn = '/run/media/ckirst/ChristophsBackuk4TB/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims';
f = io.openFile(fn);

dataset = io.readData(f, resolution=0);
dataraw = dataset[1200:1400,1200:1400,1000:1030];


# normalize data
data = dataraw.astype('float');
dmax = 0.075 * 65535;
ids = data > dmax;
data[ids] = dmax;
data /= dmax;

# plot subset od data
iplt.plotTiling(data[:,:,0:8])


# run segmentation
img = data[:,:,0:30];

import iDISCO.Segmentation.SpotDetection
import sys
self = sys.modules['iDISCO.Segmentation.SpotDetection'];



# write result










































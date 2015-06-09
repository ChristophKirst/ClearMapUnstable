# -*- coding: utf-8 -*-
"""
Processing Brain Samples - Main Example Script


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import numpy as np
import matplotlib as mpl

import iDISCO.IO.Imaris as io

from iDISCO.ImageProcessing.SpotDetection import detectCells
from iDISCO.ImageProcessing.ParallelProcessing import parallelProcessStack
from iDISCO.Visualization.Plot import plotTiling, plotOverlayLabel

from iDISCO.Utils.Timer import Timer
timer = Timer();

# open data
fn = '/home/ckirst/Science/Projects/BrainActivityMap/Data/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims';
#fn = '/run/media/ckirst/ChristophsBackuk4TB/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims';
#fn = '/home/nicolas/Windows/Nico/cfosRegistrations/Adult cfos C row 20HF 150524 - Copy.ims';
#fn = '/home/ckirst/Science/Projects/BrainActivityMap/iDISCO_2015_04/test for spots added spot.ims'



centers = parallelProcessStack(fn, chunksizemax = 40, chunksizemin = 30, chunkoverlap = 15, 
                               processes = 5, segmentation = detectCells, zrange = (0, 60));
timer.printElapsedTime("Main");


#fout = '/home/nicolas/Windows/Nico/cfosRegistrations/Adult cfos C row 20HF 150524 Points.data';
#fout = '/home/ckirst/Desktop/cfosRegistrations/Adult cfos C row 20HF 150524 Points.data';
#centers.tofile(fout);



# write result

#fout = '/home/ckirst/Science/Projects/BrainActivityMap/Data/iDISCO_2015_06/Adult cfos C row 20HF 150524 segmentation.ims';
fout = '/home/nicolas/Windows/Nico/cfosRegistrations/Adult cfos C row 20HF 150524 - Copy.ims';
h5file = io.openFile(fout);

io.writePoints(h5file, centers, mode = 'o', radius = 0.5);

h5file.close();












""" Open Data """

import iDISCO.IO.Imaris as io
import numpy as np

fn = '/home/ckirst/Science/Projects/BrainActivityMap/Data/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims';


f = io.openFile(fn);
dataset = io.readData(f, resolution=0);
print dataset.shape

dataraw = dataset[1200:1600,1200:1600,1000:1160];

print dataraw.dtype
print (dataraw.max(), dataraw.min())

datamax = np.iinfo(dataraw.dtype).max;


#plotTiling(img[:,:,0:8])
#f.close();




"""
Test between openings for speed

"""

import numpy
from mahotas import open
from scipy.ndimage.morphology import binary_opening, grey_opening
from skimage.morphology import opening, white_tophat
from skimage.filter.rank import tophat
import cv2

from iDISCO.ImageProcessing.Filter.StructureElement import structureElement
from iDISCO.Visualization.Plot import plotTiling, plotOverlayLabel

from iDISCO.Utils.Timer import Timer;


img = numpy.random.rand(2000,2000) * 65535;
img = img.astype('int')

dataraw = dataset[:,:,1160];

img = dataraw[:,:];
img.shape


t = Timer();
res = grey_opening(img, structure = structureElement('Disk', (15,15)));
t.printElapsedTime('scipy');


#t.reset();
#res2 = open(img, structureElement('Disk', (30,30)).astype('bool'));
#t.printElapsedTime('mahotas');

#t.reset();
#res2 = open(img, structureElement('Disk', (30,30)).astype('bool'));
#t.printElapsedTime('mahotas');

t.reset();
se = structureElement('Disk', (55,55)).astype('uint8');
res2 = cv2.morphologyEx(img, cv2.MORPH_OPEN, se)
t.printElapsedTime('opencv');

#from skimage.morphology import disk
#se = disk(10);

t.reset();
#res3 = opening(img.astype('uint16'), se);
res3 = white_tophat(img.astype('uint16'), se);

t.printElapsedTime('skimage');
t.reset();

x = (res3);
print (x.min(), x.max())

plotTiling(numpy.dstack((20 * img, 20 * (img - res2))))

#seems not to work correctly -> scipy gets very slow for larger images (??)




"""
Test Speed of Correlation
"""


import scipy
from iDISCO.ImageProcessing.Filter.FilterKernel import filterKernel

from skimage.feature import match_template

from iDISCO.Utils.Timer import Timer;

timer = Timer();

#DoG filter
fdog = filterKernel(ftype = 'DoG', size = [7,7,5]);
img = dataset[:,:,1000:1016];

#timer.reset();
#res = scipy.signal.correlate(img, fdog);
#timer.printElapsedTime(head = 'scipy.signal');


#res = match_template(img, fdog);
#timer.printElapsedTime(head = 'skimage');

timer.reset();
res2 = scipy.ndimage.filters.correlate(img, fdog);
timer.printElapsedTime(head = 'scipy.ndimage');
































"""

import iDISCO.IO.Imaris as io


from multiprocessing import Pool


# test parallel reading of h5py file

def readD(fnz):
    fn = fnz[0];
    z = fnz[1];
    f = io.openFile(fn, mode = 'r');
    dataset = io.readData(f, resolution=0);    
    
    img = dataset[120:130,120:200,z];
    
    f.close();
    
    return img;
    
    
fn = '/home/ckirst/Science/Projects/BrainActivityMap/iDISCO_2015_04/test for spots added spot.ims'


argdata = [(fn,x) for x in range(50)];
pool = Pool(processes = 4);

results = pool.map(readD, argdata);
print results


"""











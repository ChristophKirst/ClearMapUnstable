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

centers = parallelProcessStack(fn, chunksizemax = 40, chunksizemin = 30, chunkoverlap = 15, processes = 5, segmentation = detectCells, zrange = (0, 60));
timer.printElapsedTime("Main");


#fout = '/home/nicolas/Windows/Nico/cfosRegistrations/Adult cfos C row 20HF 150524 Points.data';
#fout = '/home/ckirst/Desktop/cfosRegistrations/Adult cfos C row 20HF 150524 Points.data';
#centers.tofile(fout);



#f = io.openFile(fn);
#dataset = io.readData(f, resolution=0);
#print dataset.shape
#f.close();
#dataraw = dataset[1200:1400,1200:1400,1000:1160];
# plot subset od data
#iplt.plotTiling(data[:,:,0:8])
#segment
#centers, imglab, imgmask = iseg.detectCells(data);

# write result

#fout = '/home/ckirst/Science/Projects/BrainActivityMap/Data/iDISCO_2015_06/Adult cfos C row 20HF 150524 segmentation.ims';
fout = '/home/nicolas/Windows/Nico/cfosRegistrations/Adult cfos C row 20HF 150524 - Copy.ims';
h5file = io.openFile(fout);

io.writePoints(h5file, centers, mode = 'o', radius = 0.5);

h5file.close();







"""
Test between mahotas and skiimage openings

"""

import numpy
from mahotas import open
from scipy.ndimage.morphology import binary_opening, grey_opening
from iDISCO.ImageProcessing.Filter.StructureElement import structureElement
from iDISCO.Visualization.Plot import plotTiling, plotOverlayLabel

from iDISCO.Utils.Timer import Timer;

t = Timer();
img = numpy.random.rand(2000,2000) * 65535;
img = img.astype('int')

res = grey_opening(img, structure = structureElement('Disk', (30,30)));

t.printElapsedTime('scipy');
t.reset();

res2 = open(img, structureElement('Disk', (30,30)).astype('bool'));
t.printElapsedTime('mahotas');

from skimage.morphology import opening

t.reset();
se = structureElement('Disk', (30,30)).astype('uint8');
res3 = opening(img.astype('uint8'), se);
t.printElapsedTime('skimage');
t.reset();

x = (res - res2);
print (x.min(), x.max())

plotTiling(numpy.dstack((img, res, res2, res3)))

#seems not to work correctly -> scipy gets very slow for larger images (??)


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











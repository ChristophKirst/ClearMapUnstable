# -*- coding: utf-8 -*-
"""
Processing Brain Samples - Main Example Script


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import numpy as np
import matplotlib as mpl

import iDISCO.Visualization.Plot as iplt
import iDISCO.Segmentation.SpotDetection as iseg
import iDISCO.IO.Imaris as io

import time

# open data
#fn = '/home/ckirst/Science/Projects/BrainActivityMap/Data/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims';
#fn = '/run/media/ckirst/ChristophsBackuk4TB/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims';
fn = '/home/nicolas/Windows/Nico/cfosRegistrations/Adult cfos C row 20HF 150524 - Copy.ims';
fn = '/home/ckirst/Science/Projects/BrainActivityMap/iDISCO_2015_04/test for spots added spot.ims'

#f = io.openFile(fn);
#dataset = io.readData(f, resolution=0);
#print dataset.shape
#f.close();
#dataraw = dataset[1200:1400,1200:1400,1000:1160];
# plot subset od data
#iplt.plotTiling(data[:,:,0:8])
#segment
#centers, imglab, imgmask = iseg.detectCells(data);

t1 = time.time();

centers = iseg.parallelProcessStack(fn, chunksizemax = 50, chunksizemin = 30, chunkoverlap = 15, processes = 5, segmentation = iseg.detectCells);

#fout = '/home/nicolas/Windows/Nico/cfosRegistrations/Adult cfos C row 20HF 150524 Points.data';
#fout = '/home/ckirst/Desktop/cfosRegistrations/Adult cfos C row 20HF 150524 Points.data';
#centers.tofile(fout);

t2 = time.time();

m, s = divmod(t2-t1, 60)
h, m = divmod(m, 60)
print "============="
print "%d:%02d:%02d" % (h, m, s)
print "============="

#f.close();


# write result

#fout = '/home/ckirst/Science/Projects/BrainActivityMap/Data/iDISCO_2015_06/Adult cfos C row 20HF 150524 segmentation.ims';
fout = '/home/nicolas/Windows/Nico/cfosRegistrations/Adult cfos C row 20HF 150524 - Copy.ims';
h5file = io.openFile(fout);

io.writePoints(h5file, centers, mode = 'o', radius = 0.5);

h5file.close();











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











# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 18:37:24 2015

@author: ckirst
"""

import os

import numpy as np
import matplotlib.pyplot as plt


import iDISCO.IO.IO as io

from iDISCO.Settings import IDISCOPath

data = io.readPoints(os.path.join(IDISCOPath, 'Data/lightsheet_line_intensty_profiles_y.csv'))


x = data[:,0]

plt.figure()
for i in range(1,data.shape[1]):
    plt.plot(x, data[:,i]);
    
    
## average
    
y = data[:,1:-1];


yn = y.copy();
yn[np.isnan(yn)] = 632;
ym = np.mean(yn, axis = 1)

plt.plot(x, ym, 'k')

## fit Guassian or r^6 poly

from scipy.optimize import curve_fit

#def f(x, a, m, s, b):
#    return a * np.exp(- (x - m)**2 / 2 / s) + b;
#
#                        
#mean = sum(ym * x)/sum(ym)
#sigma = sum(ym * (x-mean)**2)/(sum(ym))
#
#
#popt, pcov = curve_fit(f, x, ym, p0 = (1000, mean, sigma, 400))
#
#print popt
#
#def fopt(x):
#    return f(x, a = popt[0], m = popt[1], s = popt[2], b = popt[3]);
#
#plt.plot(x, map(fopt, x))


## for r^6

mean = sum(ym * x)/sum(ym)

def g(x,m,a,b,c,d):
    return a + b * (x-m)**2 + c * (x-m)**4 + d * (x-m)**6;

                        
popt, pcov = curve_fit(g, x, ym, p0 = (mean, 1, 1, 1, .1))

print popt

def gopt(x):
    return g(x, m = popt[0], a = popt[1], b = popt[2], c = popt[3], d = popt[4]);

plt.plot(x, map(gopt, x))

## create and write flatfield vector


def correctImage(fn):

    img = io.readData(fn);
    
    ds =io.dataSize(fn)
    print ds
    
    flt = map(gopt, range(0, 2560))
    #io.writePoints("/home/ckirst/Science/Projects/BrainActivityMap/Experiment/LightSheetIntensityProfile/flatfield.csv", flt)
    
    ## make flatfield
    
    fltr = flt;
    fltr.reverse();
    
    flatfield = np.zeros(ds);
    for i in range(ds[0]):
        flatfield[i,:] = fltr;
    
    import iDISCO.Visualization.Plot as dpl
    #dpl.plotTiling(np.dstack((flatfield / flatfield.max(), 2 * img.astype('float')/img.max())))
    
    
    ## correct background
    

    
    img3d = img.astype('float').copy();
    img3d.shape = (img3d.shape[0], img3d.shape[1], 1);
    
    imgc = ic.correctIllumination(img3d, flatfield);
    
    
    imgn = img.astype('float').copy();
    imgn = imgn / imgn.max();
    imgn[imgn > 0.1] = 0.1;
    
    imgcn = imgc[:,:,0].copy();
    imgcn = imgcn / imgcn.max();
    imgcn[imgcn > 0.1] = 0.1;
    
    dpl.plotTiling(np.dstack((imgn ,  imgcn)))

    io.writeData(fn[:-4] + "_c.tif", imgc)
    
    return imgc;





fn = "/home/mtllab/Desktop/vignetting/06-56-36_0_8xs3-cfos-20HFcont_UltraII_C00_xyz-Table Z1465.ome.tif"

fn = "/home/mtllab/Desktop/vignetting/13-15-16_0_8X-cfos_UltraII_C00_xyz-Table Z1367.ome.tif";
fn = "/home/mtllab/Desktop/vignetting/17-17-11_0-8xs3-cfos-20HFcont_UltraII_C00_xyz-Table Z1308.ome.tif";

fn = "/home/mtllab/Desktop/vignetting/18-29-28_0_8X-cfos_UltraII_C00_xyz-Table Z1338.ome.tif";


import iDISCO.ImageProcessing.IlluminationCorrection as ic
import iDISCO.IO.IO as io
img = io.readData(fn);
img.shape = (img.shape[0], img.shape[1], 1);
imgc = ic.correctIllumination(img);
io.writeData(fn[:-4] + "_c.tif", imgc)




#np.median(img) /  np.median(imgc)
#flatfield.sum() / flatfield.size



fns = ['/home/ckirst/Science/Projects/BrainActivityMap/Experiment/LightSheetIntensityProfile/Christoph/06-56-36_0_8xs3-cfos-20HFcont_UltraII_C00_xyz-Table Z1465.ome.tif', 
       '/home/ckirst/Science/Projects/BrainActivityMap/Experiment/LightSheetIntensityProfile/Christoph/13-15-16_0_8X-cfos_UltraII_C00_xyz-Table Z1367.ome.tif',
       '/home/ckirst/Science/Projects/BrainActivityMap/Experiment/LightSheetIntensityProfile/Christoph/17-17-11_0-8xs3-cfos-20HFcont_UltraII_C00_xyz-Table Z1308.ome.tif',
       '/home/ckirst/Science/Projects/BrainActivityMap/Experiment/LightSheetIntensityProfile/Christoph/18-29-28_0_8X-cfos_UltraII_C00_xyz-Table Z1338.ome.tif'];
       
for i in range(len(fns)):
    correctImage(fns[i])




##    

flt = map(gopt, range(0, 2560))
    
fltr = flt;
fltr.reverse();

ds = (2160, 2560);

flatfield = np.zeros(ds);
for i in range(ds[0]):
    flatfield[i,:] = fltr;

import iDISCO.Visualization.Plot as dpl
dpl.plotTiling(flatfield / flatfield.max())

io.writePoints(os.path.join(IDISCOPath, "Data/lightsheet_flatfield_correction.csv"), fltr)
    
    ## make flatfield





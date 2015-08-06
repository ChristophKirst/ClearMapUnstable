# -*- coding: utf-8 -*-
"""
Simple Interface to OME files

read data for a z-stack of tif files

Note: filename is given as regular expression

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import numpy
import matplotlib as mpl

import os
import re

import tifffile as tiff


def readFileList(filename):
    
    #get path        
    (fpath, fname) = os.path.split(filename)
    fnames = os.listdir(fpath);
    #fnames = [os.path.join(fpath, x) for x in fnames];
    
    searchRegex = re.compile(fname).search    
    fl = [ l for l in fnames for m in (searchRegex(l),) if m]  
    
    if fl == []:
        raise RuntimeError('no files found in ' + fpath + ' match ' + fname + ' !');
    
    fl.sort();
        
    return fpath, fl;


def readZRange(filename, z = all, resolution = 0):
    """Read z range from file"""

    fpath, fl = self.readFileList(filename);
    nz = len(fl);
    
    if z == all: 
        return (0,nz);
    else:
        if len(z) != 2:
            raise RuntimeError("z range has wrong format: " + str(z));
        
        if z[0] == all:
            z = (0, z[1]);
            
        if z[1] == all:
            z = (z[0], nz);
            
        if z[0] < 0 or z[0] > nz or z[1] < z[0] or z[1] > nz:
          raise RuntimeError("z range specifications out of bounds (0," + str(nz) + ") !");
           
        return z;
        
def readDataStack(filename, x = all, y = all, z = all, channel = 0, timepoint = 0, resolution = 0):
    """Read data from a single tiff stack"""
    # todo: ranges
    data = tiff.imread(filename);
    return data.transpose([1,2,0]);


def readDataFiles(filename, x = all, y = all, z = all, channel = 0, timepoint = 0, resolution = 0):
    """Read data from individual images assuming they are the z slices"""
    
    fpath, fl = self.readFileList(filename);
    nz = len(fl);
    
    #determine zrange
    if z == all:
        #assume slices go from 0 to n-1
        z = (0, nz);
    else:
        if z[0] == all:
            z = (0, z[1]);
        if z[1] == all:
            z = (z[0], nz);
        nz = z[1] - z[0];
    
    fn = os.path.join(fpath, fl[0]);    
    img = mpl.pyplot.imread(fn);
    
    if x == all:
        nx = img.shape[0];
        x = (0, nx);
    else:
        if x[0] == all:
            x = (0, x[1]);
        if x[1] == all:
            x = (x[0], img.shape[0]);    
        nx = x[1] - x[0];
            
    if y == all:
        ny = img.shape[1];
        y = (0, ny);
    else:
        if y[0] == all:
            y = (0, y[1]);
        if y[1] == all:
            y = (y[0], img.shape[1]);        
        ny = y[1] - y[0];        
    
    img = numpy.zeros((nx,ny,nz), dtype = img.dtype);
    
    for iz in range(z[0], z[1]):
        fn = os.path.join(fpath, fl[iz]);  
        imgraw = mpl.pyplot.imread(fn);
        img[:,:,iz - z[0]] = imgraw[x[0]:x[1], y[0]:y[1]];
    
    return img;

     
def readData(filename, x = all, y = all, z = all, channel = 0, timepoint = 0, resolution = 0):
    """Read a tiff stack as multiple form single or multiple images"""
    
    if os.path.exists(filename):
         readDataFunction = self.readDataStack;
    else:
         readDataFunction = self.readDataFiles;
    
    return readDataFunction(filename, x = x, y = y, z = z, channel = channel, resolution = resolution);



def writeData(img, fileheader):
    nz = img.shape[2];
    #vmin = numpy.iinfo(img.dtype).min
    #vmax = numpy.iinfo(img.dtype).max
    for i in range(nz):
        #mpl.pyplot.imsave(fileheader + "%.4d" % i + ".tif", img[:,:,i], vmin = vmin, vmax = vmax, cmap = mpl.pyplot.cm.gray);
        tiff.imsave(fileheader + "%.4d" % i + ".tif", img[:,:,i])


def writeDataStack(filename, image):
    """Writes volumetric data to a single tiff stack assuming (y,x,z) array representation"""

    tiff.imsave(filename, image.transpose([2,0,1]));
 




def test():    
    """Test OME module"""  
    fn = os.path.split(self.__file__);
    fn = os.path.join(fn[0], '../Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z\d{4}.ome.tif');
    
    print "Searching images of form: " + fn;
    
    fp, fl = self.readFileList(fn);
    print "Found " + str(len(fl)) + " images!"
    
    zr = self.readZRange(fn);
    print "Found z range: " + str(zr);
    
    print "Loading full z-stack..."
    img = self.readData(fn);  
    print "Image size: " + str(img.shape)
    
    print "Loading subset of z-stack..."
    img = self.readData(fn, z = (1,3), x = (500,600), y = (600,all));  
    print "Image size: " + str(img.shape)
        

if __name__ == "__main__":
    
    self.test();
    

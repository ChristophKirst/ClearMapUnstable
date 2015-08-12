# -*- coding: utf-8 -*-
"""
Simple Interface to stack saved as a list of files

Note: filename is given as regular expression

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import numpy
import tifffile as tiff

import iDISCO.IO.IO as io


def dataSize(filename, **args):
    t = tiff.TiffFile(filename);
    d3 = len(t.pages);
    d2 = t.pages[0].shape;
    d2 = (d2[0], d2[1]);
    if d3 > 1:
        dims = d2 + (d3,);
    else:
        dims =  d2;
    
    return io.dataSizeFromDataRange(dims, **args);

def dataZSize(filename, z = all, **args):
    t = tiff.TiffFile(filename);
    d3 = len(t.pages);
    if d3 > 1:
        return io.toDataSize(d3, r = z);
    else:
        return None;


def readData(filename, x = all, y = all, z = all, **args):
    """Read data from a single tif image or stack"""
    
    dsize = self.dataSize(filename);
    if len(dsize) == 2:
        data = tiff.imread(filename, key = 0);
        return io.dataToRange(data.transpose([0,1]), x = x, y = y);
    else:
        if z is all:
            data = tiff.imread(filename);
            return io.dataToRange(data.transpose([1,2,0]), x = x, y = y, z = all);
        else: #optimize for z ranges
            dsize = io.dataSizeFromDataRange(dsize, x = x, y = y, z = z);
            t = tiff.Tifffile(filename);
            p = t.pages[0];
            data = numpy.zeros(dsize, dtype = p.dtype);
            rz = io.toDataRange(dsize[2], r = z);
            for i in range(rz[0], rz[1]):
                xydata = t.pages[i].asarray();
                data[:,:,i-rz[0]] = io.dataToRange(xydata, x = x, y = y);
            
        return data


def writeData(filename, data):
    d = len(data.shape);
    
    if d == 2:
        tiff.imsave(filename, data);
    elif d == 3:   
        tiff.imsave(filename, data.transpose([2,0,1]));
    else:
        raise RuntimeError('writing multiple channel data to tif not supported');
    

def test():    
    """Test TIF module"""  
    #import iDISCO.IO.TIF as self
    
    from iDISCO.Parameter import iDISCOPath
    import os
    import numpy

    basedir = iDISCOPath();
    fn = os.path.join(basedir,'Test/Data/Tif/test.tif')  

    data = numpy.random.rand(20,50,10);
    data[5:15, 20:45, 2:9] = 0;
    data = 20 * data;
    data = data.astype('int32');

    #reload(self)
    print "writing raw image to: " + fn;    
    self.writeData(fn, data);

    print "Loading raw image from: " + fn;
    img = self.readData(fn);  
    print "Image size: " + str(img.shape)
    
    diff = img - data;
    print (diff.max(), diff.min())

    
    fn = os.path.join(basedir,'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z1000.ome.tif')        
    img = self.readData(fn);  
    print "Image size: " + str(img.shape)
    
    #dataSize
    print "dataSize  is %s" % str(self.dataSize(fn))
    print "dataZSize is %s" % str(self.dataZSize(fn))
    
    print "dataSize  is %s" % str(self.dataSize(fn, x = (10,20)))
    print "dataZSize is %s" % str(self.dataZSize(fn))
        

if __name__ == "__main__":
    
    self.test();
    

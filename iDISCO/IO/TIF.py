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
    #d2 = (d2[0], d2[1]);
    d2 = (d2[1], d2[0]);
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
    #print "dsize %s" % str(dsize);    
    
    if len(dsize) == 2:
        data = tiff.imread(filename, key = 0);
        #print "data.shape %s" % str(data.shape);        
        
        return io.dataToRange(data.transpose([1,0]), x = x, y = y);
        #return io.dataToRange(data, x = x, y = y);
        
    else:
        if z is all:
            data = tiff.imread(filename);
            if data.ndim == 2:
                # data = data
                data = data.transpose([1,0]);
            elif data.ndim == 3:
                #data = data.transpose([1,2,0]);
                data = data.transpose([2,1,0]);
            elif data.ndim == 4: # multi channel image
                #data = data.transpose([1,2,0,3]);
                data = data.transpose([2,1,0,3]);
            else:
                raise RuntimeError('readData: dimension %d not supproted!' % data.ndim)
            
            return io.dataToRange(data, x = x, y = y, z = all);
        
        else: #optimize for z ranges
            ds = io.dataSizeFromDataRange(dsize, x = x, y = y, z = z);
            t = tiff.TiffFile(filename);
            p = t.pages[0];
            data = numpy.zeros(ds, dtype = p.dtype);
            rz = io.toDataRange(dsize[2], r = z);
            
            print "test"
            print rz;
            print dsize            
            
            for i in range(rz[0], rz[1]):
                xydata = t.pages[i].asarray();
                #data[:,:,i-rz[0]] = io.dataToRange(xydata, x = x, y = y);
                data[:,:,i-rz[0]] = io.dataToRange(xydata.transpose([1,0]), x = x, y = y);
            
            return data


def writeData(filename, data):
    d = len(data.shape);
    
    if d == 2:
        #tiff.imsave(filename, data);
        tiff.imsave(filename, data.transpose([1,0]));
    elif d == 3:   
        #tiff.imsave(filename, data.transpose([2,0,1]));
        tiff.imsave(filename, data.transpose([2,1,0]));
    elif d == 4:
        #tiffile (z,y,x,c)
        t = tiff.TiffWriter(filename, bigtiff = True);
        #t.save(data.transpose([2,0,1,3]), photometric = 'minisblack',  planarconfig = 'contig');
        t.save(data.transpose([2,1,0,3]), photometric = 'minisblack',  planarconfig = 'contig')
        t.close();    
    else:
        raise RuntimeError('writing multiple channel data to tif not supported');
    
    return filename;
    

def copyData(source, sink):
    io.copyFile(source, sink);


def test():    
    """Test TIF module"""  
    #import iDISCO.IO.TIF as self
    
    from iDISCO.Settings import iDISCOPath
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
    
    
    print "Loading raw image from %s with limited z range: " % fn;
    img = self.readData(fn, z = (3,8));  
    print "Image size: " + str(img.shape)
    
    diff = img - data[:,:,3:8];
    print (diff.max(), diff.min())

    
    fn = os.path.join(basedir,'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z1000.ome.tif')        
    
    ds = self.dataSize(fn);
    print "Image size form dataSiZe: " + str(ds)    
    
    
    img = self.readData(fn);  
    print "Image size: " + str(img.shape)
    
    #dataSize
    print "dataSize  is %s" % str(self.dataSize(fn))
    print "dataZSize is %s" % str(self.dataZSize(fn))
    
    print "dataSize  is %s" % str(self.dataSize(fn, y = (10,20)))
    print "dataZSize is %s" % str(self.dataZSize(fn))
        
    #test writing multi channel image
    x = numpy.random.rand(50,100,30,4) * 10;
    x = x.astype('float32');
    self.writeData(os.path.join(basedir,'Test/Data/Tif/composite.tif'), x)

    y = self.readData(os.path.join(basedir,'Test/Data/Tif/composite.tif'));
    y.shape
        
    diff = x - y;
    print (diff.max(), diff.min())


if __name__ == "__main__":
    
    self.test();
    



# Test various tif libs:
#import numpy;
#import iDISCO.IO.IO as io
#x = numpy.random.rand(30,50,100,4) * 10;
#x = x.astype('float32');
#
#io.writeData('test.tif', x.astype('float32'))
#
#import tifffile as tiff
#
#tiff.imsave('test.tif', d)   
#
#t = tiff.TiffWriter('test.tif', bigtiff = True);
#
#t.save(x.astype('int16'), photometric = 'minisblack',  planarconfig = 'contig')
##t.save(x.astype('int16'), photometric = 'minisblack',  planarconfig = None)
#t.close();
#
#d = tiff.imread('Composite.tif')
#
#from libtiff import TIFF
#tif = TIFF.open('Composite.tif', mode='r')
#image = tif.read_image()
#image.shape
#
#from PIL import Image
#im = Image.open('Composite.tif');
#img = numpy.array(im);
#img.shape
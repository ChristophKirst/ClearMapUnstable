# -*- coding: utf-8 -*-
"""
Simple Interface to image stack saved as a list of files

Note: filename is given as regular expression

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import numpy

import os
import re

import iDISCO.IO.IO as io


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


def dataSize(filename, **args):
    fp, fl = self.readFileList(filename);
    nz = len(fl);
    
    d2 = io.dataSize(os.path.join(fp, fl[0]));
    if not len(d2) == 2:
        raise RuntimeError("FileList: importing multiple files of dim %d not supported!" % len(d2));
    
    dims = d2 + (nz,);
    return io.dataSizeFromDataRange(dims, **args);
   
   
def dataZSize(filename, z = all, **args):
    fp, fl = self.readFileList(filename);
    nz = len(fl);
    return io.toDataSize(nz, r = z);



def readDataFiles(filename, x = all, y = all, z = all, **args):
    """Read data from individual images assuming they are the z slices"""   
    
    fpath, fl = self.readFileList(filename);
    nz = len(fl);
    
    #read first image to get data size and type
    rz = io.toDataRange(nz, r = z);
    sz = io.toDataSize(nz, r = z);
    fn = os.path.join(fpath, fl[rz[0]]); 
    img = io.readData(fn, x = x, y = y);
    nxy = img.shape;
    data = numpy.zeros(nxy + (sz,), dtype = img.dtype);
    data[:,:,0] = img;


    for i in range(rz[0]+1, rz[1]):
        fn = os.path.join(fpath, fl[i]);    
        data[:,:,i-rz[0]] = io.readData(fn, x = x, y = y);
    
    return data;

     
def readData(filename, **args):
    """Read image stack from single or multiple images"""
    
    if os.path.exists(filename):
         return io.readData(filename, **args);
    else:
         return self.readDataFiles(filename, **args);


def splitFileExpression(filename):
    
    searchRegex = re.compile('.*\\\\d\{(?P<digit>\d)\}.*').search
    m = searchRegex(filename); 
    
    if not m is None:
        digits = int(m.group('digit'));
        searchRegex = re.compile('.*(?P<replace>\\\\d\{\d\}).*').search
        m = searchRegex(filename);       
        fs = filename.split(m.group('replace'));
        if not len(fs) == 2:
            raise RuntimeError('FileList: more than a single replacement for z indexing!');
        fileheader = fs[0];
        fileext    = fs[1];
        
    else: #fileheader is given
        digits = 4;
        fileheader = filename;
        fileext    = '.tif';
        
    digitfrmt = "%." + str(digits) + "d";        
        
    return (fileheader, fileext, digitfrmt);



def writeData(filename, data):
    """Write image stack to single or multiple image files"""
    #check for the \d{xx} part of the regular expression -> if not assume file header

    (fileheader, fileext, digitfrmt) = self.splitFileExpression(filename);

    d = len(data.shape);
    if d == 2:
        fname = fileheader + (digitfrmt % 0) + fileext;
        io.writeData(fname, data);
        return fname;
    else:
        nz = data.shape[2];
        for i in range(nz):
            fname = fileheader + (digitfrmt % i) + fileext;
            io.writeData(fname, data[:,:,i]);
        return filename;


def copyData(source, sink):
    
    (fileheader, fileext, digitfrmt) = self.splitFileExpression(sink);
    
    fp, fl = self.readFileList(source);
    
    for i in range(len(fl)):
        io.copyFile(os.path.join(fp, fl[i]), fileheader + (digitfrmt % i) + fileext);
    
    return sink


def test():    
    """Test FileList module"""  
    import iDISCO.IO.FileList as self
    
    from iDISCO.Parameter import iDISCOPath
    import os
    import numpy

    basedir = iDISCOPath();
    fn = os.path.join(basedir,'Test/Data/FileList/test\d{4}.tif')  

    data = numpy.random.rand(20,50,10);
    data[5:15, 20:45, 2:9] = 0;
    data = 20 * data;
    data = data.astype('int32');

    reload(self)
    print "writing raw image to: " + fn;    
    self.writeData(fn, data);

    print "Loading raw image from: " + fn;
    img = self.readData(fn);  
    print "Image size: " + str(img.shape)
    
    diff = img - data;
    print (diff.max(), diff.min())

    
    fn = os.path.join(basedir,'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z\d{4}.ome.tif')        
    
    fp, fl = self.readFileList(fn);
    print "Found " + str(len(fl)) + " images!"    
    
    #dataSize
    print "dataSize  is %s" % str(self.dataSize(fn))
    print "dataZSize is %s" % str(self.dataZSize(fn))
    
    print "dataSize  is %s" % str(self.dataSize(fn, x = (10,20)))
    print "dataZSize is %s" % str(self.dataZSize(fn))
    
    
    img = self.readData(fn, z = (17,all));  
    print "Image size: " + str(img.shape)    
    
    import iDISCO.Visualization.Plot as plt
    plt.plotTiling(img)
        

if __name__ == "__main__":
    
    self.test();
    

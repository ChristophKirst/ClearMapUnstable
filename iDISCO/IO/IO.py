# -*- coding: utf-8 -*-
"""
IO interface to read microscope and point data

Notes:
    - Numerical arrays represent images as (y,x,z) or (y,x) for 2d images
    - Center of cells / points are represented in the same format (y,x,z) or (y,x)
    - Supported image file extensions are given in the imageFileExtensions list
    - Supported point files extensions are given in the pointFileExtensions list

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import os
import re
import numpy
import importlib

pointFileExtensions = ["csv", "txt", "vtk", "ims"];

pointFileTypes = ["CSV", "VTK", "Imaris"];

pointFileExtensionToType = {"csv":"CSV", "txt":"CSV", "vtk":"VTK", "ims":"Imaris"};

imageFileExtensions = ["tif", "tiff", "mdh", "raw", "ims", "nrrd"]

imageFileTypes = ["FileList", "TIF", "RAW", "NRRD", "Imaris"]

imageFileExtensionToType = { "tif" : "TIF", "tiff" : "TIF", "raw" : "RAW", "mhd" : "RAW", "nnrd": "NRRD", "ims" : "Imaris"}


##############################################################################
# Basic file querries
##############################################################################

def fileExtension(filename):
    """Returns file extension if exists"""
    
    fext = filename.split('.');
    if len(fext) < 2:
        return None;
    else:
        return fext[-1];
        
    
def isFile(filename):
    """Checks if filename is a real file, returns false if it is directory or regular expression"""
  
    if os.path.exists(filename):
        if os.path.isdir(filename):
            return False;
        else:
            return True;
    else:
        return False;
    

def isFileExpression(filename):
    if self.isFile(filename):
        return False;
    else:
        searchRegex = re.compile('.*\\\\d\{(?P<digit>\d)\}.*').search
        m = searchRegex(filename);
        if m is None:
            return False;
        else:
            return True;
            

def pointFileNameToType(filename):
    """Finds type of a point file"""
    
    #if not self.isFile(filename):
    #    raise RuntimeError("Cannot find point file %s" % filename);
    #else:
    fext = self.fileExtension(filename);
    if fext in pointFileExtensions:
        return pointFileExtensionToType[fext];
    else:
       raise RuntimeError("Cannot read point file %s with extension %s" % (filename, fext));     



def imageFileNameToType(filename):
    """Finds type of an image file"""
    
    if self.isFileExpression(filename):
        return "FileList";
    else:
        fext = self.fileExtension(filename);
        if fext in  imageFileExtensions:
            return imageFileExtensionToType[fext];
        else:
           raise RuntimeError("Cannot read image file %s with extension %s" % (filename, fext));


def imageFileNameToModule(filename):
    """Return the module that handles io for the image file filename"""
    ft = self.imageFileNameToType(filename);
    return importlib.import_module("iDISCO.IO." + ft);


def pointFileNameToModule(filename):
    """Return the module that handles io for the point file filename"""
    ft = self.pointFileNameToType(filename);
    return importlib.import_module("iDISCO.IO." + ft);



##############################################################################
# Read / Write Data
##############################################################################

    
def dataSize(filename, x = all, y = all, z = all, **args):
    """Returns size of the array needed to read the data in filename"""        
    mod = self.imageFileNameToModule(filename);
    return mod.dataSize(filename, x = x, y = y, z = z, **args);
    
def dataZSize(filename, z = all, **args):
    """Returns size of the array in the third dimension, None if 2D data"""
    mod = self.imageFileNameToModule(filename);
    return mod.dataZSize(filename, z = z, **args);

def toDataRange(size, r = all):
    if r is all:
        return (0,size);
    else:
        if r[0] is all:
            r = (0, r[1]);
        if r[0] < 0:
            if -r[0] > size:
                r = (0, r[1]);
            else:
                r = (size + r[0], r[1]);
        if r[0] > size:
            r = (size, r[1]);
            
        if r[1] is all:
            r = (r[0], size);
        if r[1] < 0:
            if -r[1] > size:
                r = (r[0], 0);
            else:
                r = (r[0], size + r[1]);
        if r[1] > size:
            r = (r[0], size);
        
        if r[0] > r[1]:
            r = (r[0], r[0]);
        
        return r;
        
def toDataSize(size, r = all):
    dr = self.toDataRange(size, r = r);
    return int(dr[1] - dr[0]);

def dataSizeFromDataRange(dsize, x = all, y = all, z = all, **args):
    dsize = list(dsize);
    n = len(dsize);
    if n > 0:
        dsize[0] = self.toDataSize(dsize[0], r = y);
    if n > 1:
        dsize[1] = self.toDataSize(dsize[1], r = x);
    if n > 2:
        dsize[2] = self.toDataSize(dsize[2], r = z);
    
    return tuple(dsize);

def dataToRange(data, x = all, y = all, z = all, **args):
    dsize = data.shape;
    n = len(dsize);
    rr = [];
    if n > 0:
        rr.append(self.toDataRange(dsize[0], r = y));
    if n > 1:
        rr.append(self.toDataRange(dsize[1], r = x));
    if n > 2:
        rr.append(self.toDataRange(dsize[2], r = z));
    
    if n == 1:
        return data[rr[0][0] : rr[0][1]];
    elif n == 2:
        return data[rr[0][0] : rr[0][1], rr[1][0] : rr[1][1]];
    elif n == 3:
        return data[rr[0][0] : rr[0][1], rr[1][0] : rr[1][1], rr[2][0] : rr[2][1]];
    else:
        raise RuntimeError('dataToRangt dealing with images larger than 4d not supported');
            


def readData(filename, **args):
    """Read data from one of the supported formats"""
    
    if filename is None:
        return None;   
    elif isinstance(filename, basestring):
        mod = self.imageFileNameToModule(filename);
        return mod.readData(filename, **args);
    elif isinstance(filename, numpy.ndarray ):
        return filename;
    else:
        raise RuntimeError('readData: cannot infer format of the requested data/file.');

def writeData(filename, data, **args):
    """Write data to file filename"""
    
    if filename is None:
        return None;
    
    mod = self.imageFileNameToModule(filename);    
    return mod.writeData(filename, data, **args);



##############################################################################
# Read / Write Points
##############################################################################

   
def readPoints(filename, **args):
    """Read a list of points from csv or vtk"""
    if filename == None:
        return
    
    mod = self.pointFilenameToModule(filename);
    return mod.readPoints(filename, **args);

   
def writePoints(filename, points, **args):
    """Write a list of points to csv, vtk or ims files"""
    if filename == None:
        return
    
    mod = self.pointFilenameToModule(filename);
    return mod.writePoints(filename, points, **args);


 
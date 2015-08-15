# -*- coding: utf-8 -*-
"""
IO interface to read microscope and point data

Notes:
    - Numerical arrays represent data as (y,x,z) or (y,x) for 2d images
    - Center of cells / points are represented in the same format (y,x,z) or (y,x)
    - Supported image file extensions are given in the imageFileExtensions list
    - Supported point files extensions are given in the pointFileExtensions list
    - all routines accesing data or data properties accept filenames and numpy arrays

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import os
import re
import numpy
import importlib
import shutil

pointFileExtensions = ["csv", "txt", "vtk", "ims"];

pointFileTypes = ["CSV", "VTK", "Imaris"];

pointFileExtensionToType = {"csv" : "CSV", "txt" : "CSV", "vtk" : "VTK", "ims" : "Imaris"};

dataFileExtensions = ["tif", "tiff", "mdh", "raw", "ims", "nrrd"]

dataFileTypes = ["FileList", "TIF", "RAW", "NRRD", "Imaris"]

dataFileExtensionToType = { "tif" : "TIF", "tiff" : "TIF", "raw" : "RAW", "mhd" : "RAW", "nrrd": "NRRD", "ims" : "Imaris"}


##############################################################################
# Basic file querries
##############################################################################

def fileExtension(source):
    """Returns file extension if exists"""
    
    if not isinstance(source, basestring):
        return None;
    
    fext = source.split('.');
    if len(fext) < 2:
        return None;
    else:
        return fext[-1];
        
    
def isFile(source):
    """Checks if filename is a real file, returns false if it is directory or regular expression"""
    
    if not isinstance(source, basestring):
        return False;
  
    if os.path.exists(source):
        if os.path.isdir(source):
            return False;
        else:
            return True;
    else:
        return False;
    

def isFileExpression(source):
    """Checks if filename is a regular expression denoting a file list"""    
    
    if not isinstance(source, basestring):
        return False;
    
    if self.isFile(source):
        return False;
    else:
        searchRegex = re.compile('.*\\\\d\{(?P<digit>\d)\}.*').search
        m = searchRegex(source);
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
       raise RuntimeError("Cannot determine type of point file %s with extension %s" % (filename, fext));     


def dataFileNameToType(filename):
    """Finds type of an data file"""
    
    if self.isFileExpression(filename):
        return "FileList";
    else:
        fext = self.fileExtension(filename);
        if fext in dataFileExtensions:
            return dataFileExtensionToType[fext];
        else:
           raise RuntimeError("Cannot determine type of data file %s with extension %s" % (filename, fext));


def dataFileNameToModule(filename):
    """Return the module that handles io for the data file filename"""
    ft = self.dataFileNameToType(filename);
    return importlib.import_module("iDISCO.IO." + ft);


def pointFileNameToModule(filename):
    """Return the module that handles io for the point file filename"""
    ft = self.pointFileNameToType(filename);
    return importlib.import_module("iDISCO.IO." + ft);



##############################################################################
# Read / Write Data
##############################################################################

    
def dataSize(source, x = all, y = all, z = all, **args):
    """Returns size of the array needed to read the data in filename"""  
    
    if isinstance(source, basestring):
        mod = self.dataFileNameToModule(source);
        return mod.dataSize(source, x = x, y = y, z = z, **args);
    elif isinstance(source, numpy.ndarray):
        return source.shape;
    else:
        raise RuntimeError("dataSize: argument not a string or array!");
            
    
def dataZSize(source, z = all, **args):
    """Returns size of the array in the third dimension, None if 2D data"""
      
    if isinstance(source, basestring):
        mod = self.dataFileNameToModule(source);
        return mod.dataZSize(source, z = z, **args);
    elif isinstance(source, numpy.ndarray):
        if len(source.shape) > 2: 
            return self.toDataSize(source.shape[2], r = z);
        else:
            return None;   
    else:
        raise RuntimeError("dataZSize: argument not a string or array!");        


def toDataRange(size, r = all):
    """Converts range r to numeric range (min,max) given the full array size"""
    if r is all:
        return (0,size);
    
    if isinstance(r, int) or isinstance(r, float):
        r = (r, r +1);
      
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
    """Converts full size to actual size given range r"""
    dr = self.toDataRange(size, r = r);
    return int(dr[1] - dr[0]);


def dataSizeFromDataRange(dsize, x = all, y = all, z = all, **args):
    """Converts full data size to actual size given ranges for x,y,z"""
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
    d = len(dsize);
    rr = [];
    if d > 0:
        rr.append(self.toDataRange(dsize[0], r = y));
    if d > 1:
        rr.append(self.toDataRange(dsize[1], r = x));
    if d > 2:
        rr.append(self.toDataRange(dsize[2], r = z));
    if d > 4:
        raise RuntimeError('pointsToRange: dimension %d to big' % d);
    
    if d == 1:
        return data[rr[0][0] : rr[0][1]];
    elif d == 2:
        return data[rr[0][0] : rr[0][1], rr[1][0] : rr[1][1]];
    elif d == 3:
        return data[rr[0][0] : rr[0][1], rr[1][0] : rr[1][1], rr[2][0] : rr[2][1]];
    elif d == 4:
        return data[rr[0][0] : rr[0][1], rr[1][0] : rr[1][1], rr[2][0] : rr[2][1], :];   
    else:
        raise RuntimeError('dataToRange dealing with images of dimension %d not supported' % d);
            

def readData(source, **args):
    """Read data from one of the supported formats or numpy array"""
    
    if source is None:
        return None;   
    elif isinstance(source, basestring):
        mod = self.dataFileNameToModule(source);
        return mod.readData(source, **args);
    elif isinstance(source, numpy.ndarray ):
        return self.dataToRange(source, **args);
    else:
        raise RuntimeError('readData: cannot infer format of the requested data/file.');

def writeData(sink, data, **args):
    """Write data to file filename"""
    
    if sink is None: #dont write to disk but return the data
        return data;
    
    mod = self.dataFileNameToModule(sink);    
    return mod.writeData(sink, data, **args);


def copyFile(source, sink):
    shutil.copy(source, sink);

def copyData(source, sink):
    mod = self.dataFileNameToModule(source);
    return mod.copyData(source, sink);


def transformData(source, sink, **args):
    """Transforms data in source to sink"""
    
    if source is None:
        return None;   
    
    elif isinstance(source, basestring):
        if sink is None:        
            return self.readData(source, **args);
        elif isinstance(sink, basestring):
            if args == {} and self.dataFileNameToType(source) == self.dataFileNameToType(sink):
                return self.copyData(source, sink);
            else:
                data = self.readData(source, **args);
                return self.writeData(sink, data);
        else:
            raise RuntimeError('transformData: unknown sink!');
            
    elif isinstance(source, numpy.ndarray):
        if sink is None:
            return self.dataToRange(source, **args);
        elif isinstance(sink,  basestring):
            data = self.dataToRange(source, **args);
            return self.writeData(sink, data);
        else:
            raise RuntimeError('transformData: unknown sink!');
 
 
def toMultiChannelData(*args):
    """Concatenate single channel arrays to one multi channel array"""
    data = numpy.array(args);
    return data.rollaxis(data, 0, data.ndim);


##############################################################################
# Read / Write Points
##############################################################################

def points(points):
    if isinstance(points, tuple):
        return points[0];
    else:
        return points;
        
def intensities(points):
    if isinstance(points, tuple):
        return points[1];
    else:
        return None;


def pointShiftFromRange(dataSize, x = all, y = all, z = all, **args):
    """Calculate shift of points given a specific range restriction"""
    
    if isinstance(dataSize, basestring):
        dataSize = self.dataSize(dataSize);
    
    d = len(dataSize.shape);
    rr = [];
    if d > 0:
        rr.append(self.toDataRange(dataSize[0], r = y));
    if d > 1:
        rr.append(self.toDataRange(dataSize[1], r = x));
    if d > 2:
        rr.append(self.toDataRange(dataSize[2], r = z));
    if d > 3 or d < 1:
        raise RuntimeError('shiftFromRange: dimension %d to big' % d);
    
    return [r[0] for r in rr];


def pointsToRange(points, dataSize = all, x = all, y = all, z = all, shift = False, **args):
    """Restrict points to a specific range"""
        
    d = points.shape[1];
    
    if dataSize is all:
        dataSize = points.max(axis=1);
    elif isinstance(dataSize, basestring):
        dataSize = self.dataSize(dataSize);
    
    rr = [];
    if d > 0:
        rr.append(self.toDataRange(dataSize[0], r = y));
    if d > 1:
        rr.append(self.toDataRange(dataSize[1], r = x));
    if d > 2:
        rr.append(self.toDataRange(dataSize[2], r = z));
    if d > 3 or d < 1:
        raise RuntimeError('pointsToRange: dimension %d to big' % d);
    
    #ids = numpy.zeros(n, dtype = 'bool');
    if d > 0:
        ids = numpy.logical_and(points[:,0] >= rr[0][0], points[:,0] <= rr[0][1]);
    if d > 1:
        ids = numpy.logical_and(numpy.logical_and(ids, points[:,1] >= rr[1][0]), points[:,1] <= rr[1][1]);
    if d > 2:
        ids = numpy.logical_and(numpy.logical_and(ids, points[:,2] >= rr[2][0]), points[:,2] <= rr[2][1]);
        
    points = points[ids, :];
    
    if shift:
        sh = [r[0] for r in rr];
        points = points - sh;
    
    return points;


def readPoints(source, **args):
    """Read a list of points from csv or vtk"""
    if source is None:
        return None;
    elif isinstance(source, basestring):
        mod = self.pointFileNameToModule(source);
        return mod.readPoints(source, **args);
    elif isinstance(source, numpy.ndarray):
        return self.pointsToRange(source, **args);
    else:
        raise RuntimeError('readPoints: cannot infer format of the requested data/file.');


def writePoints(sink, points, **args):
    """Write a list of points to csv, vtk or ims files"""
    if sink is None:
        return points;
    
    mod = self.pointFileNameToModule(sink);
    return mod.writePoints(sink, points, **args);


 
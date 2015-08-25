# -*- coding: utf-8 -*-
"""
IO interface to read microscope and point data

Notes:
    - Numerical arrays represent data as (x,z,z) or (x,y) for 2d images
    - Center of cells / points are represented in the same format (x,y,z) or (x,y)
    - Supported image file extensions are given in the imageFileExtensions list
    - Supported point files extensions are given in the pointFileExtensions list
    - all routines accesing data or data properties accept filenames and numpy arrays
    
    - Note many image reading routines in other packages read images as (y,x,z) or (y,x) arrays

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

pointFileExtensions = ["csv", "txt", "npy", "vtk", "ims"];

pointFileTypes = ["CSV", "NPY", "VTK", "Imaris"];

pointFileExtensionToType = {"csv" : "CSV", "txt" : "CSV", "npy" : "NPY", "vtk" : "VTK", "ims" : "Imaris"};

dataFileExtensions = ["tif", "tiff", "mhd", "raw", "ims", "nrrd"]

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
            
def createDirectory(filename):
    dirname, fname = os.path.split(filename);
    if not os.path.exists(dirname):
        os.makedirs(dirname);
    

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
        dsize[0] = self.toDataSize(dsize[0], r = x);
    if n > 1:
        dsize[1] = self.toDataSize(dsize[1], r = y);
    if n > 2:
        dsize[2] = self.toDataSize(dsize[2], r = z);
    
    return tuple(dsize);

def dataToRange(data, x = all, y = all, z = all, **args):
    dsize = data.shape;
    d = len(dsize);
    rr = [];
    if d > 0:
        rr.append(self.toDataRange(dsize[0], r = x));
    if d > 1:
        rr.append(self.toDataRange(dsize[1], r = y));
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


def convertData(source, sink, **args):
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



def sagittalToCoronalData(source, sink = None):
    """Change from saggital to coronal orientation"""
      
    source = self.readData(source);
    d = source.ndim;
    if d < 3:
        raise RuntimeError('sagittalToCoronalData: 3d image requited!');
    
    tp = range(d);
    tp[0:3] = [2,0,1];
    source = source.transpose(tp);
    source = source[::-1];
    #source = source[::-1,:,:];
    return self.writeData(sink, source);



##############################################################################
# Read / Write Points
##############################################################################


def pointsToCoordinates(points):
    if isinstance(points, tuple):
        return points[0];
    else:
        return points;
        
def pointsToProperties(points):
    if isinstance(points, tuple) and len(points) > 1:
        return points[1];
    else:
        return None;
        
def pointsToCoordinatesAndProperties(points):
    if isinstance(points, tuple):
        if len(points) == 0:
            return (None, None);
        elif len(points) == 1:
            return (points[0], None);
        elif len(points) == 2:
            return points;
        else:
            raise RuntimeError('points not a tuple of 0 to 2 elements!');
    else:
        return (points, None);


def pointsToCoordinatesAndPropertiesFileNames(filename, propertiesPostfix = '_intensities', **args):
    if isinstance(filename, basestring):
        return (filename, filename[:-4] + propertiesPostfix + filename[-4:])
    elif isinstance(filename, tuple):
        if len(filename) == 1:
            if filename[0] is None:
                return (None, None);
            elif isinstance(filename[0], basestring):
                return (filename[0], filename[0][:-4] + propertiesPostfix + filename[0][-4:]);
            else:
                raise RuntimeError('pointsFilenames: invalid filename specification!');
        elif len(filename) == 2:
            return filename;
        else:
            raise RuntimeError('pointsFilenames: invalid filename specification!');
    elif filename is None:
        return (None, None)
    else:
        raise RuntimeError('pointsFilenames: invalid filename specification!');




def pointShiftFromRange(dataSize, x = all, y = all, z = all, **args):
    """Calculate shift of points given a specific range restriction"""
    
    if isinstance(dataSize, basestring):
        dataSize = self.dataSize(dataSize);
    dataSize = list(dataSize);
    
    d = len(dataSize);
    rr = [];
    if d > 0:
        rr.append(self.toDataRange(dataSize[0], r = x));
    if d > 1:
        rr.append(self.toDataRange(dataSize[1], r = y));
    if d > 2:
        rr.append(self.toDataRange(dataSize[2], r = z));
    if d > 3 or d < 1:
        raise RuntimeError('shiftFromRange: dimension %d to big' % d);
    
    return [r[0] for r in rr];

def pointsToRange(points, dataSize = all, x = all, y = all, z = all, shift = False,  **args):
    """Restrict points to a specific range"""

    if x is all and y is all and z is all:
        return points;

    istuple = isinstance(points, tuple);    
    (points, properties) = pointsToCoordinatesAndProperties(points);
    
    if points is None:
        if istuple:
            return (points, properties);
        else:
            return points;

    if not isinstance(points, numpy.ndarray):
        raise RuntimeError('pointsToRange: points not None or numpy array!');
        
        
    d = points.shape[1];
    
    if dataSize is all:
        dataSize = points.max(axis=0);
    elif isinstance(dataSize, basestring):
        dataSize = self.dataSize(dataSize);
    
    rr = [];
    if d > 0:
        rr.append(self.toDataRange(dataSize[0], r = x));
    if d > 1:
        rr.append(self.toDataRange(dataSize[1], r = y));
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
    
    if not properties is None:
        properties = properties[ids];
    
    if istuple:
        return (points, properties);
    else:
        return points;


def readPoints(source, **args):
    """Read a list of points from csv or vtk"""
    
    istuple = isinstance(source, tuple);

    if source is None:
        source = (None, None);
    elif isinstance(source, numpy.ndarray):
        source = (source, None);
    elif isinstance(source, basestring):
        source = (source, None);
    elif isinstance(source, tuple):
        if len(source) == 0:
            source = (None, None);
        elif len(source) == 1: 
            if source[0] is None:
                source = (None, None);
            elif isinstance(source[0], numpy.ndarray):
                source = (source[0], None);
            elif isinstance(source[0], basestring):
                source = pointsToCoordinatesAndPropertiesFileNames(source, **args);
            else:
                raise RuntimeError('readPoints: cannot infer format of the requested data/file.');       
        elif len(source) == 2:
            if not((source[0] is None or isinstance(source[0], basestring) or isinstance(source[0], numpy.ndarray)) and 
                   (source[1] is None or isinstance(source[1], basestring) or isinstance(source[0], numpy.ndarray))):
               raise RuntimeError('readPoints: cannot infer format of the requested data/file.');
        else:
            raise RuntimeError('readPoints: cannot infer format of the requested data/file.');
    else:
        raise RuntimeError('readPoints: cannot infer format of the requested data/file.'); 
           
    if source[0] is None:
        points = None;
    elif isinstance(source[0], numpy.ndarray):
        points = source[0];
    elif isinstance(source[0], basestring):
        mod = self.pointFileNameToModule(source[0]);
        points = mod.readPoints(source[0]);

    if source[1] is None:
        properties = None;
    elif isinstance(source[1], numpy.ndarray):
        properties = source[1];
    elif isinstance(source[1], basestring):
        mod = self.pointFileNameToModule(source[1]);
        properties = mod.readPoints(source[1]);
        
    if istuple:
        return self.pointsToRange((points, properties), **args);
    else:
        return self.pointsToRange(points, **args);


def writePoints(sink, points, **args):
    """Write a list of points to csv, vtk or ims files"""
    
    
    #todo: make clean independent of return of two results -> io.wrtiePoints -> take care of pairs: (points,intensities)
    istuple = isinstance(sink, tuple);    
    
    if sink is None:
        sink = (None, None);
    elif isinstance(sink, basestring):
        sink = (sink, None);
    elif isinstance(sink, tuple):
        if len(sink) == 0:
            sink = (None, None);
        elif len(sink) == 1:
            if sink[0] is None:
                sink = (None, None);
            elif isinstance(sink, basestring):
                sink = pointsToCoordinatesAndPropertiesFileNames(sink, **args);
            else:
                raise RuntimeWarning('sink not well defined!')
        elif len(sink) == 2:
            if not((sink[0] is None or isinstance(sink[0], basestring)) and (sink[1] is None or isinstance(sink[1], basestring))):
                raise RuntimeWarning('sink not well defined!')
        else:
            raise RuntimeWarning('sink not well defined!')
    else:
        raise RuntimeWarning('sink not well defined!')

    
    (points, properties) = pointsToCoordinatesAndProperties(points); 
    if sink[0] is None:
        retpoints = points;
    else:
        mod = self.pointFileNameToModule(sink[0]);
        retpoints = mod.writePoints(sink[0], points);
    
    if sink[1] is None:
        retproperties = properties;
    else:
        mod = self.pointFileNameToModule(sink[1]);
        retproperties = mod.writePoints(sink[1], properties);
        
    if istuple:
        return (retpoints, retproperties);
    else:
        return retpoints;
        

 
 
 
 
 
 
 
 
 





# -*- coding: utf-8 -*-
"""
IO interface to read microscope data

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

from iDISCO.IO import CSV, Imaris, OME, VTK, RAW

import os

def readZRange(filename, z = all, resolution = 0):
    """Determine z range of file filename"""
        
    #get extension
    fext = filename.split('.');
    if len(fext) < 2:
        raise RuntimeError('readZRange: cannot infer file format without extension');
        
    fext = fext[-1];
    
    if fext == "ims":
        readZRangeFunction = Imaris.readZRange;
    elif fext == "tif":
        readZRangeFunction = OME.readZRange;
    else:
        raise RuntimeError('readZRange: file format ' + fext + ' not ims or tif');

    return readZRangeFunction(filename, z = z, resolution = resolution);


def readData(filename, x = all, y = all, z = all, channel = 0, timepoint = 0, resolution = 0):
    """Read data from one of the supported formats"""
    
    #get extension
    fext = filename.split('.');
    if len(fext) < 2:
        raise RuntimeError('readData: cannot infer file format without extension');
        
    fext = fext[-1];
    
    if fext == "ims":
        readDataFunction = Imaris.readData;
    elif fext == "tif":
        if os.path.exists(filename):
            readDataFunction = OME.readDataStack;
        else:
            readDataFunction = OME.readData;
    elif fext == "mhd":
        readDataFunction = RAW.readData;
    else:
        raise RuntimeError('readData: file format ' + fext + ' not supported');

    return readDataFunction(filename, x = x, y = y, z = z, channel = channel, resolution = resolution);
    
def redDataStack(filename):
    """Read data from a single tiff stack"""
    return OME.readDataStack(filename);
    

def writePoints(filename, points):
    """Write a list of points to csv, vtk or ims files"""
    if filename == None:
        return
    
    #get extension
    fext = filename.split('.');
    if len(fext) < 2:
        raise RuntimeError('readData: cannot infer file format without extension');
        
    fext = fext[-1];
    
    if fext == "ims":
        writePointsFunction = Imaris.writePoints;
        points = Imaris.transformToImaris(points);
    elif fext == "csv":
        writePointsFunction = CSV.writePoints;
    elif fext == "vtk":
        writePointsFunction = VTK.writePoints;
    else:
        raise RuntimeError('readData: file format ' + fext + ' not ims, vtk or csv');

    return writePointsFunction(filename, points);
    
    
def readPoints(filename):
    """Read a list of points from csv or vtk"""
    if filename == None:
        return
    
    #get extension
    fext = filename.split('.');
    if len(fext) < 2:
        raise RuntimeError('readData: cannot infer file format without extension');
        
    fext = fext[-1];
    
    if fext == "ims":
        readPointsFunction = CSV.readPoints
    elif fext == "csv":
        readPointsFunction = CSV.readPoints;
    elif fext == "vtk":
        readPointsFunction = VTK.readPoints;
    else:
        raise RuntimeError('readData: file format ' + fext + ' not ims, vtk or csv');

    return readPointsFunction(filename);

    
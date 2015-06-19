# -*- coding: utf-8 -*-
"""
IO interface to read microscope data

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import iDISCO.IO as io

def readZRange(filename, z = all, resolution = 4):
    """Determine z range of file filename"""
        
    #get extension
    fext = filename.split('.');
    if len(fext) < 2:
        raise RuntimeError('readZRange: cannot infer file format without extension');
        
    fext = fext[-1];
    
    if fext == "ims":
        readZRangeFunction = io.Imaris.readZRange;
    elif fext == "tif":
        readZRangeFunction = io.OME.readZRange;
    else:
        raise RuntimeError('readZRange: file format ' + fext + ' not ims or tif');

    return readZRangeFunction(filename, z = z, resolution = resolution);

    


def readData(filename, x = all, y = all, z = all, channel = 0, timepoint = 0, resolution = 4):
    """Read data from one of the supported formats"""
    
    #get extension
    fext = filename.split('.');
    if len(fext) < 2:
        raise RuntimeError('readData: cannot infer file format without extension');
        
    fext = fext[-1];
    
    if fext == "ims":
        readDataFunction = io.Imaris.readData;
    elif fext == "tif":
        readDataFunction = io.OME.readData;
    else:
        raise RuntimeError('readData: file format ' + fext + ' not ims or tif');

    return readDataFunction(filename, x = x, y = y, z = z, channel = channel, resolution = resolution);
    
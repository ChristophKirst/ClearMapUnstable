# -*- coding: utf-8 -*-
"""
Simple Interface to read RAW/MHD files e.g. created by elastix

Created on Thu Jun  4 14:37:06 2015

Todo: read subsets efficiently

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import os
import numpy

import vtk
from vtk.util.numpy_support import vtk_to_numpy

import iDISCO.IO.IO as io


def dataSize(filename, **args):
    """Read data size from raw/mhd image"""
    
    imr = vtk.vtkMetaImageReader()
    imr.SetFileName(filename);
    imr.Update()

    im = imr.GetOutput()
    dims = im.GetDimensions();
    dims = list(dims);
    dims[0:2] = [dims[1], dims[0]];
    dims = tuple(dims);
    return io.dataSizeFromDataRange(dims, **args);

    
def dataZSize(filename, z = all, **args):
    """Read data size from raw/mhd image"""
    
    imr = vtk.vtkMetaImageReader()
    imr.SetFileName(filename);
    imr.Update()

    im = imr.GetOutput()
    dims = im.GetDimensions()
    
    if len(dims) > 2:
        return io.toDataSize(dims[2], r = z);
    else:
        return None;

    

def readData(filename, x = all, y = all, z = all):
    """Read data from raw/mhd image"""
    
    imr = vtk.vtkMetaImageReader()
    imr.SetFileName(filename);
    imr.Update()

    im = imr.GetOutput()
    dims = im.GetDimensions()
    sc = im.GetPointData().GetScalars()
    img = vtk_to_numpy(sc)
    #print img.shape

    dims = (dims[2], dims[1], dims[0]);
    img = img.reshape(dims)
    img = img.transpose([1,2,0]);
    
    return io.dataToRange(img, x = x, y = y, z = z);


def writeHeader(filename, meta_dict):
    header = ''
    # do not use tags = meta_dict.keys() because the order of tags matters
    tags = ['ObjectType','NDims','BinaryData',
       'BinaryDataByteOrderMSB','CompressedData','CompressedDataSize',
       'TransformMatrix','Offset','CenterOfRotation',
       'AnatomicalOrientation',
       'ElementSpacing',
       'DimSize',
       'ElementType',
       'ElementDataFile',
       'Comment','SeriesDescription','AcquisitionDate','AcquisitionTime','StudyDate','StudyTime']
    for tag in tags:
        if tag in meta_dict.keys():
            header += '%s = %s\n'%(tag,meta_dict[tag])
    f = open(filename,'w')
    f.write(header)
    f.close()


def writeRawData(filename, data):
    """ Write the data into a raw format file. Big endian is always used. """    
    rawfile = open(filename,'wb');
    d = len(data.shape);
    if d <= 2:
        data.tofile(rawfile);
    elif d == 3:
        data.transpose([2,0,1]).tofile(rawfile);
    elif d== 4:
        data.transpose([3,2,0,1]).tofile(rawfile);
    else:
        raise RuntimeError('writeRawData: image dimension %d not supported!' % d);
    
    rawfile.close();
    
    return filename;


def writeData(filename, data, **args):
    
    fext = io.fileExtension(filename);
    if fext == "raw":
        fname = filename[:-3] + 'mhd';
    else:
        fname = filename;
    
    assert(fname[-4:]=='.mhd')
    
    meta_dict = {}
    meta_dict['ObjectType'] = 'Image'
    meta_dict['BinaryData'] = 'True'
    meta_dict['BinaryDataByteOrderMSB'] = 'False'

    numpy_to_datatype   = {numpy.dtype('int8')    : "MET_CHAR",
                           numpy.dtype('uint8')   : "MET_UCHAR",
                           numpy.dtype('int16')   : "MET_SHORT", 
                           numpy.dtype('uint16')  : "MET_USHORT",
                           numpy.dtype('int32')   : "MET_INT",
                           numpy.dtype('uint32')  : "MET_UINT",
                           numpy.dtype('int64')   : "MET_LONG",
                           numpy.dtype('uint64')  : "MET_ULONG",
                           numpy.dtype('float32') : "MET_FLOAT", 
                           numpy.dtype('float64') : "MET_DOUBLE",
                           }    
                           
    dtype = data.dtype;    
    meta_dict['ElementType'] = numpy_to_datatype[dtype];
    
    #fix arrays represented as (y,x,z)
    dsize = list(data.shape);    
    dsize[0:2] = [dsize[1],dsize[0]];
    meta_dict['NDims'] = str(len(dsize))
    meta_dict['DimSize'] = ' '.join([str(i) for i in dsize])
    meta_dict['ElementDataFile'] = os.path.split(fname)[1].replace('.mhd','.raw')
    self.writeHeader(fname, meta_dict)

    pwd = os.path.split(fname)[0]
    if pwd:
        data_file = pwd +'/' + meta_dict['ElementDataFile']
    else:
        data_file = meta_dict['ElementDataFile']

    self.writeRawData(data_file, data)
    
    return fname;


def copyData(source, sink):
    sourceExt = io.fileExtension(source);
    sinkExt   = io.fileExtension(sink);
    
    sources = [source]; 
    sinks = [];
    
    if sourceExt == 'raw':
        sources.append(source[:-3] + 'mhd');
        
        if sinkExt == 'raw':
            sinks.append(sink);
            sinks.append(sink[:-3] + 'mhd');
        elif sinkExt == 'mhd':
            sinks.append(sink[:-3] + 'raw');
            sinks.append(sink);
        else:
            raise RuntimeError('copyData: sink extension %s not raw or mhd' % sinkExt);
    
    elif sourceExt == 'mhd':
        sources.append(source[:-3] + 'raw');
        
        if sinkExt == 'raw':
            sinks.append(sink[:-3] + 'mhd');
            sinks.append(sink);
        elif sinkExt == 'mhd':
            sinks.append(sink);
            sinks.append(sink[:-3] + 'raw');
        else:
            raise RuntimeError('copyData: sink extension %s not raw or mhd' % sinkExt);
        
    for i in range(2):
        io.copyData(sources[i], sinks[i]);
    
    return sink;



def test():    
    from iDISCO.Parameter import iDISCOPath
    import os
    import numpy
    
    """Test RAW module"""
    basedir = iDISCOPath();
    fn = os.path.join(basedir, 'Test/Data/Raw/test.mhd')

    data = numpy.random.rand(20,50,10);
    data[5:15, 20:45, 2:9] = 0;

    #reload(self)
    print "writing raw image to: " + fn;    
    self.writeData(fn, data);

    print "Loading raw image from: " + fn;
    img = self.readData(fn);  
    print "Image size: " + str(img.shape)
    
    diff = img - data;
    print (diff.max(), diff.min())

    #some uint type
    print "writing raw image to: " + fn;    
    udata = data * 10;
    udata = udata.astype('uint16');
    self.writeData(fn, udata);

    print "Loading raw image from: " + fn;
    img = self.readData(fn);  
    print "Image size: " + str(img.shape)
    
    diff = img - udata;
    print (diff.max(), diff.min())
    
    
    #dataSize
    print "dataSize  is %s" % str(self.dataSize(fn))
    print "dataZSize is %s" % str(self.dataZSize(fn))
    

if __name__ == "__main__":
    self.test();
    


#VTK / ITK versions:

#def writeData(filename, data, **args):
#    
#    vtkarray = numpy_to_vtk(data);
#    #assert isinstance(v, vtk.vtkImageData)
#    writer = vtk.vtkMetaImageWriter();
#    
#    fext = io.fileExtension(filename);
#    if fext == "raw":
#        fname = filename[:-3] + 'mhd';
#    else:
#        fname = filename;
#        
#    writer.SetFileName(fname)
#    writer.SetInput(vtkarray)
#    writer.Write()


#   from medpy.io import save    
#   save(data, filename);   
#   return filename;
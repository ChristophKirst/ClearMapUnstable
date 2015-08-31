# -*- coding: utf-8 -*-
"""
Resample data to match Allen Brain Atlas

Notes:
 - uses iDISCO convetion of (x,y) array representation of an image
 - points are assumed to be given in iDISCO x,y,z representation
 - parameters such as resolutions and orientation or datasize are assumed to be given as x,y,z as when viewed e.g. in imagej 
 - a permutation is defined as permuting the axes and then inverting the axes 

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import multiprocessing  
import tempfile
#import matplotlib.pyplot as plt
import os
import math
import numpy
import shutil
import cv2

from iDISCO.Utils.ProcessWriter import ProcessWriter;
import iDISCO.IO.IO as io
import iDISCO.IO.FileList as fl


def fixOrientation(orientation):
    """Substitue Orientation Strings"""
    
    if orientation is None:
        return None;
        
    #fix named representations
    if orientation == 'Left':
        orientation = (1,2,3);
    if orientation == 'Right':
        orientation = (-1,2,3);    
    
    return orientation;


def inverseOrientation(orientation):
    """Returns the inverse permuation of the permutation per taking axis inversions into account"""
    if orientation is None:
        return None;
    
    n = len(orientation);
    iper = list(orientation);
    
    #permutation is defined as permuting the axes and then axis inversion
    for i in range(n):
        if orientation[i] < 0:
            iper[int(abs(orientation[i])-1)] = -(i + 1);
        else:
            iper[int(abs(orientation[i])-1)] = (i + 1);
    
    return tuple(iper)


def orientationToPermuation(orientation):
    orientation = fixOrientation(orientation);
    if orientation is None:
        return (0,1,2);
    else:
        return tuple(int(abs(i))-1 for i in orientation);


def orientResolution(resolution, orientation):
    if resolution is None:
        return None;
    
    per = orientationToPermuation(orientation);
    #print orientation, per, resolution
    return tuple(resolution[i] for i in per);
    

def orientResolutionInverse(resolution, orientation):
    if resolution is None:
        return None;
    
    per = orientationToPermuation(inverseOrientation(orientation));
    return tuple(resolution[i] for i in per);

 
def orientDataSize(dataSize, orientation):
    return orientResolution(dataSize, orientation);
 
def orientDataSizeInverse(dataSize, orientation):
    return orientResolutionInverse(dataSize, orientation); 
 
 
def resampleDataSize(dataSizeSource, dataSizeSink = None, resolutionSource = None, resolutionSink = None, orientation = None):
    """Calculate scaling factors and data sizes for resampling"""

    orientation = fixOrientation(orientation);    
    
    #determine data sizes if not specified
    if dataSizeSink is None:
        if resolutionSource is None or resolutionSink is None:
            raise RuntimeError('resampleDataSize: data size and resolutions not defined!');
        
        #orient resolution of source to resolution of sink to get sink data size
        resolutionSourceO = orientResolution(resolutionSource, orientation);
        
        #calculate scaling factor
        dataSizeSink = tuple([int(math.ceil(dataSizeSource[i] *  resolutionSourceO[i]/resolutionSink[i])) for i in range(len(dataSizeSource))]);        
        
    #print dataSizeSink, "ds sink"
    
    if dataSizeSource is None:
        if resolutionSource is None or resolutionSink is None:
            raise RuntimeError('resampleDataSize: data size and resolutions not defined!');
        
        #orient resolution of source to resolution of sink to get sink data size
        resolutionSourceO = orientResolution(resolutionSource, orientation);
        
        #calculate source data size
        dataSizeSource = tuple([int(math.ceil(dataSizeSink[i] *  resolutionSink[i]/resolutionSource[i])) for i in range(len(dataSizeSink))]);  
        dataSizeSource = orientDataSizeInverse(dataSizeSource);
        
    #print dataSizeSource, "ds source"
        
    #calculate effecive resolutions
    if resolutionSource is None:
        if resolutionSink is None:
            resolutionSource = (1,1,1);
        else:
            dataSizeSourceO = orientDataSize(dataSizeSource, orientation);
            resolutionSource = tuple(float(dataSizeSink[i]) / dataSizeSourceO[i] * resolutionSink[i] for i in range(len(dataSizeSource)));
            resolutionSource = orientResolutionInverse(resolutionSource, orientation);
    
    #print resolutionSource, "res source sink"
    
    dataSizeSourceO = orientDataSize(dataSizeSource, orientation);
    
    
    resolutionSourceO = orientResolution(resolutionSource, orientation);
    resolutionSink = tuple(float(dataSizeSourceO[i]) / float(dataSizeSink[i]) * resolutionSource[i] for i in range(len(dataSizeSource)));
    
    #print dataSizeSource, dataSizeSink, resolutionSource, resolutionSink 
    
    return dataSizeSource, dataSizeSink, resolutionSource, resolutionSink  




def fixInterpolation(interpolation):
    if interpolation == 'nn' or interpolation is None or interpolation == cv2.INTER_NEAREST:
        interpolation = cv2.INTER_NEAREST;
    else:
        interpolation = cv2.INTER_LINEAR;
        
    return interpolation;
        


def resampleXY(source, dataSizeSink, sink = None, interpolation = 'linear', out = sys.stdout):
    """Resample a 2d image"""   
    
    #out.write("Input: %s Output: " % (inputFile, soutputFile))
    data = io.readData(source);
    dataSize = data.shape;
    
    #print dataSize, dataSizeSink    
    
    if data.ndim != 2:
        raise RuntimeError('resampleXY: expects 2d image source, found %dd' % data.ndim)
    #print sagittalImageSize;
    
    #dataSizeSink = tuple([int(math.ceil(dataSize[i] *  resolutionSource[i]/resolutionSink[i])) for i in range(2)]);
    out.write(("resampleData: Imagesize: %d, %d " % (dataSize[0], dataSize[1])) + ("Resampled Imagesize: %d, %d" % (dataSizeSink[0], dataSizeSink[1])))
    #out.write(("resampleData: Imagesize: %d, %d " % dataSize) + ("Resampled Imagesize: %d, %d" % (outputSize[1], outputSize[0])))
    
    # note: cv2.resize reverses x-Y axes
    interpolation = fixInterpolation(interpolation)
    sinkData = cv2.resize(data,  (dataSizeSink[1], dataSizeSink[0]), interpolation = interpolation);
    #sinkData = cv2.resize(data,  outputSize);
    #sinkData = scipy.misc.imresize(sagittalImage, outputImageSize, interp = 'bilinear'); #normalizes images -> not usefull for stacks !
    
    #out.write("resampleData: resized Image size: %d, %d " % sinkData.shape)
    
    return io.writeData(sink, sinkData);


def resampleXYParallel(arg):
    """Resampling function to use for parallel resmapling of mage slices"""
    
    fileSource = arg[0];
    fileSink = arg[1];
    dataSizeSink = arg[2];
    interpolation = arg[3];
    ii = arg[4];
    nn = arg[5];
    
    pw = ProcessWriter(ii);
    pw.write("resampleData: resampling in XY: image %d / %d" % (ii, nn))
    
    data = numpy.squeeze(io.readData(fileSource, z = ii));
    self.resampleXY(data, sink = fileSink, dataSizeSink = dataSizeSink, interpolation = interpolation, out = pw);




def resampleData(source, sink = None,  orientation = None, dataSizeSink = None, resolutionSource = (4.0625, 4.0625, 3), resolutionSink = (25, 25, 25), 
                 processingDirectory = None, processes = 1, cleanup = True, verbose = True, interpolation = 'linear', **args):
    """Resample data in source in resolution and orientation"""   
    """Notes: resolutions are assumed to be given for the axes of the intrinsic orientation of the data and reference as when viewed by matplotlib or ImageJ,
              orientation: permuation of 1,2,3 with potential sign, indicating which axes map onto the reference axes, 
              negative sign indicates reversal of that particular axes        
    """
        
    orientation = fixOrientation(orientation);
    
    if isinstance(dataSizeSink, basestring):
        dataSizeSink = io.dataSize(dataSizeSink);

    #orient actual resolutions onto reference resolution    
    dataSizeSource = io.dataSize(source);
        
    dataSizeSource, dataSizeSink, resolutionSource, resolutionSink = resampleDataSize(dataSizeSource = dataSizeSource, dataSizeSink = dataSizeSink, 
                                                                                      resolutionSource = resolutionSource, resolutionSink = resolutionSink, orientation = orientation);
    
    dataSizeSinkI = orientDataSizeInverse(dataSizeSink, orientation);
    
    #print dataSizeSource, dataSizeSink, resolutionSource, resolutionSink, dataSizeSinkI
    
     
    #rescale in x y in parallel
    if processingDirectory == None:
        processingDirectory = tempfile.mkdtemp();     
        
    interpolation = fixInterpolation(interpolation);
     
    nZ = dataSizeSource[2];
    pool = multiprocessing.Pool(processes=processes);
    argdata = [];
    for i in range(nZ):
        argdata.append( (source, os.path.join(processingDirectory, 'resample_%04d.tif' % i), dataSizeSinkI, interpolation, i, nZ) );  
        #print argdata[i]
    pool.map(self.resampleXYParallel, argdata);
    
    #rescale in z
    fn = os.path.join(processingDirectory, 'resample_%04d.tif' % 0);
    data = io.readData(fn);
    zImage = numpy.zeros((dataSizeSinkI[0], dataSizeSinkI[1], nZ), dtype = data.dtype);    
    for i in range(nZ):
        if verbose and i % 10 == 0:
            print "resampleData; reading %d/%d" % (i, nZ);
        fn = os.path.join(processingDirectory, 'resample_%04d.tif' % i);
        zImage[:,:, i] = io.readData(fn);

    
    resampledData = numpy.zeros(dataSizeSinkI, dtype = zImage.dtype);

    for i in range(dataSizeSinkI[0]):
        if verbose and i % 25 == 0:
            print "resampleData: processing %d/%d" % (i, dataSizeSinkI[0])
        #resampledImage[:, iImage ,:] =  scipy.misc.imresize(zImage[:,iImage,:], [resizedZAxisSize, sagittalImageSize[1]] , interp = 'bilinear'); 
        #cv2.resize takes reverse order of sizes !
        resampledData[i ,:, :] =  cv2.resize(zImage[i,:,:], (dataSizeSinkI[2], dataSizeSinkI[1]), interpolation = interpolation);
        #resampledData[i ,:, :] =  cv2.resize(zImage[i,:, :], (dataSize[1], resizedZSize));
    

    #account for using (z,y,x) array representation -> (y,x,z)
    #resampledData = resampledData.transpose([1,2,0]);
    #resampledData = resampledData.transpose([2,1,0]);
    
    if cleanup:
        shutil.rmtree(processingDirectory);

    if not orientation is None:
        
        #reorient
        per = orientationToPermuation(orientation);
        resampledData = resampledData.transpose(per);
    
        #reverse orientation after permuting e.g. (-2,1) brings axis 2 to first axis and we can reorder there
        if orientation[0] < 0:
            resampledData = resampledData[::-1, :, :];
        if orientation[1] < 0:
            resampledData = resampledData[:, ::-1, :]; 
        if orientation[2] < 0:
            resampledData = resampledData[:, :, ::-1];
        
        #bring back from y,x,z to z,y,x
        #resampledImage = resampledImage.transpose([2,0,1]);
    if verbose:
        print "resampleData: resampled data size: " + str(resampledData.shape)  
    
    if sink == []:
        if io.isFileExpression(source):
            sink = os.path.split(source);
            sink = os.path.join(sink[0], 'resample_\d{4}.tif');
        elif isinstance(source, basestring):
            sink = source + '_resample.tif';
        else:
            raise RuntimeError('resampleData: automatic sink naming not supported for non string source!');
    
    return io.writeData(sink, resampledData);
    
    
    


def resampleDataInverse(sink, source = None, dataSizeSource = None, orientation = None, resolutionSource = (4.0625, 4.0625, 3), resolutionSink = (25, 25, 25), 
                        processingDirectory = None, processes = 1, cleanup = True, verbose = True, interpolation = 'linear', **args):
    """Resample data inversely to resampleData routine"""
    
    #orientation
    orientation = fixOrientation(orientation);
    
    #assume we can read data fully into memory
    resampledData = io.readData(sink);

    dataSizeSink = resampledData.shape;
    
    if isinstance(dataSizeSource, basestring):
        dataSizeSource = io.dataSize(dataSizeSource);

    dataSizeSource, dataSizeSink, resolutionSource, resolutionSink = resampleDataSize(dataSizeSource = dataSizeSource, dataSizeSink = dataSizeSink, 
                                                                                      resolutionSource = resolutionSource, resolutionSink = resolutionSink, orientation = orientation);

    #print (dataSizeSource, dataSizeSink, resolutionSource, resolutionSink )
    
    dataSizeSinkI = orientDataSizeInverse(dataSizeSink, orientation);
    
    
    #flip axes back and permute inversely
    if not orientation is None:
        if orientation[0] < 0:
            resampledData = resampledData[::-1, :, :];
        if orientation[1] < 0:
            resampledData = resampledData[:, ::-1, :]; 
        if orientation[2] < 0:
            resampledData = resampledData[:, :, ::-1];

        
        #reorient
        peri = self.inverseOrientation(orientation);
        peri = orientationToPermuation(peri);
        resampledData = resampledData.transpose(peri);
    
    # upscale in z
    interpolation = fixInterpolation(interpolation);
    
    resampledDataXY = numpy.zeros((dataSizeSinkI[0], dataSizeSinkI[1], dataSizeSource[2]), dtype = resampledData.dtype);    
    
    for i in range(dataSizeSinkI[0]):
        if verbose and i % 25 == 0:
            print "resampleDataInverse: processing %d/%d" % (i, dataSizeSinkI[0])

        #cv2.resize takes reverse order of sizes !
        resampledDataXY[i ,:, :] =  cv2.resize(resampledData[i,:,:], (dataSizeSource[2], dataSizeSinkI[1]), interpolation = interpolation);

    # upscale x, y in parallel
    
    if io.isFileExpression(source):
        files = source;
    else:
        if processingDirectory == None:
            processingDirectory = tempfile.mkdtemp();   
        files = os.path.join(sink[0], 'resample_\d{4}.tif');
    
    io.writeData(files, resampledDataXY);
    
    nZ = dataSizeSource[2];
    pool = multiprocessing.Pool(processes=processes);
    argdata = [];
    for i in range(nZ):
        argdata.append( (source, fl.fileExperssionToFileName(files, i), dataSizeSource, interpolation, i, nZ) );  
    pool.map(self.resampleXYParallel, argdata);
    
    if io.isFileExpression(source):
        return source;
    else:
        data = io.convertData(files, source);
        
        if cleanup:
            shutil.rmtree(processingDirectory);
        
        return data;
    



def resamplePoints(pointSource, pointSink = None, dataSizeSource = None, dataSizeSink = None, orientation = None, resolutionSource = (4.0625, 4.0625, 3), resolutionSink = (25, 25, 25), **args):
    """Resample Points to map from original data to the coordinates of the resampled image"""
    
    #fix (y,x,z) image array representation
    #resolutionSource, resolutionSink = self.fixResolutions(resolutionSource, resolutionSink);
    
    orientation = self.fixOrientation(orientation);

    #datasize of data source
    if isinstance(dataSizeSource, basestring):
        dataSizeSource = io.dataSize(dataSizeSource);
    
    dataSizeSource, dataSizeSink, resolutionSource, resolutionSink = resampleDataSize(dataSizeSource = dataSizeSource, dataSizeSink = dataSizeSink, 
                                                                                      resolutionSource = resolutionSource, resolutionSink = resolutionSink, orientation = orientation);

    points = io.readPoints(pointSource);

    dataSizeSinkI = orientDataSizeInverse(dataSizeSink, orientation);
    #resolutionSinkI = orientResolutionInverse(resolutionSink, orientation);
        
    #scaling factors
    scale = [float(dataSizeSource[i]) / float(dataSizeSinkI[i]) for i in range(3)];
    #print scale
    
    repoints = points.copy();
    for i in range(3):    
        repoints[:,i] = repoints[:,i] / scale[i];
               
    #permute for non trivial orientation
    if not orientation is None:
        per = orientationToPermuation(orientation);
        repoints = repoints[:,per];
        
        for i in range(3):
            if orientation[i] < 0:
                repoints[:,i] = dataSizeSink[i] - repoints[:,i];
      
    return io.writePoints(pointSink, repoints);


     
def resamplePointsInverse(pointSource, pointSink = None, dataSizeSource = None, dataSizeSink = None, orientation = None, resolutionSource = (4.0625, 4.0625, 3), resolutionSink = (25, 25, 25), **args):
    """Resample Points to map from the coordinates of the resampled image to the original data"""
       
    orientation = fixOrientation(orientation);
    
    #datasize of data source
    if isinstance(dataSizeSource, basestring):
        dataSizeSource = io.dataSize(dataSizeSource);
    
    dataSizeSource, dataSizeSink, resolutionSource, resolutionSink = resampleDataSize(dataSizeSource = dataSizeSource, dataSizeSink = dataSizeSink, 
                                                                                      resolutionSource = resolutionSource, resolutionSink = resolutionSink, orientation = orientation);
            
    points = io.readPoints(pointSource);
    
    dataSizeSinkI = orientDataSizeInverse(dataSizeSink, orientation);
    #resolutionSinkI = orientResolutionInverse(resolutionSink, orientation);
        
    #scaling factors
    scale = [float(dataSizeSource[i]) / float(dataSizeSinkI[i]) for i in range(3)];
    #print scale

    rpoints = points.copy();    
    
    #invert axis inversion and permutations    
    if not orientation is None:
        #invert permuation
        iorientation = self.inverseOrientation(orientation);
        per = orientationToPermuation(iorientation);
        rpoints = rpoints[:,per];
        
        for i in range(3):
            if iorientation[i] < 0:
                rpoints[:,i] = dataSizeSink[i] - rpoints[:,i];
    
    #scale points
    for i in range(3):   
        rpoints[:,i] = rpoints[:,i] * scale[i];    
    
    return io.writePoints(pointSink, rpoints);




def test():
    """Test Resampling Module"""
    import iDISCO.Alignment.Resampling as self
    reload(self)
    from iDISCO.Settings import IDISCOPath 
    import iDISCO.IO.IO as io
    import os, numpy

    fn = os.path.join(IDISCOPath, 'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z\d{4}.ome.tif');
    outfn = os.path.join(IDISCOPath, "Test/Data/Resampling/test.mhd")
    
    print "Making resampled stack " + outfn
    print "source datasize %s" % str(io.dataSize(fn));
    data = self.resampleData(fn, sink = None, resolutionSource = (1,1,1), orientation = (1,2,3), resolutionSink = (10,10,2));
    print data.shape
    io.writeData(outfn, data)   

    data = self.resampleData(fn, sink = None, dataSizeSink = (50,70,10), orientation = (1,2,3));
    print data.shape
    io.writeData(outfn, data)   


    dataSizeSource, dataSizeSink, resolutionSource, resolutionSink = self.resampleDataSize(dataSizeSource = (100,200, 303), dataSizeSink = None, 
                                                                                      resolutionSource = (1,1,1), resolutionSink = (5,5,5), orientation = (1,2,3));

    print dataSizeSource, dataSizeSink, resolutionSource, resolutionSink
    

    points = numpy.array([[0,0,0], [1,1,1], io.dataSize(fn)]);
    points = points.astype('float')
    pr = self.resamplePoints(points, dataSizeSource = fn, dataSizeSink = (50,70,10), orientation = (1,2,3))
    print pr

    pri = self.resamplePointsInverse(pr, dataSizeSource = fn, dataSizeSink = (50,70,10), orientation = (-1,2,3))
    print pri


    result = self.resampleDataInverse(outfn, os.path.join(IDISCOPath, 'Test/Data/OME/resample_\d{4}.ome.tif'), dataSizeSource = fn);

if __name__ == "__main__":
    self.test();
    




#    
#def dataSize(imageFilePattern):
#    """Determine full size from raw data in (x,y,z) order (not the array format (y,x,z))"""
#
#    if os.path.exists(imageFilePattern): # single file
#        tf = tiff.TiffFile(imageFilePattern);
#        shape = tf.series[0]['shape'];
#        return (shape[1], shape[2], shape[0])
#    
#    imageDirectory, listOfImages = self.readFileList(imageFilePattern);
#    nz = len(listOfImages);
#    
#    if nz == 0:
#        raise RuntimeError("dataSize: no files match: %s" % imageFilePattern);
#    
#    imagefile = os.path.join(imageDirectory, listOfImages[0]);
#    sagittalImage = plt.imread(imagefile);# reads as y-x
#    
#    return  (sagittalImage.shape[1], sagittalImage.shape[0], nz)
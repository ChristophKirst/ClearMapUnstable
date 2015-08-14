# -*- coding: utf-8 -*-
"""
Resample data to match Allen Brain Atlas

Notes:
 - uses iDISCO convetion of (y,x) array representation of an image
 - points are assumed to be given in iDISCO y,x,z representation
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


def fixOrientation(orientation):
    """Fix orientation as images are represented as (y,x,z) arrays but user provides orientation as premutation of x,y,z"""
    
    if orientation is None:
        return None;
        
    #fix named representations
    if orientation == 'Left':
        orientation = (1,2,3);
    if orientation == 'Right':
        orientation = (-1,2,3);
            
    #fix (y,x) image representation / take into account sign too
    o = list(orientation);
    for i in range(3):
        if int(abs(orientation[i])) == 1:
            pos1 = i;
        if int(abs(orientation[i])) == 2:
            pos2 = i;
    
    if orientation[pos2] < 0:
        o[pos2] = -abs(orientation[pos1]);
    else:
        o[pos2] = abs(orientation[pos1]);
        
    if orientation[pos1] < 0:
        o[pos1] = -abs(orientation[pos2]);
    else:
        o[pos1] = abs(orientation[pos2]);
        
    orientation = (o[1],o[0],o[2]); 
    
    return orientation;


def fixResolutions(resolutionSource, resolutionSink):
    """Fix resolutions as images are represented as (y,x,z) arrays but the user provides resolutions in x,y,z"""

    rd = (resolutionSource[1], resolutionSource[0], resolutionSource[2]);
    rr = (resolutionSink[1], resolutionSink[0], resolutionSink[2]);
    
    return (rd,rr);  
   
def fixDataSize(dataSize):
    """Fix data sizes as images are represented as (y,x,z) arrays but the user provides datasize as x,y,z"""
    return (dataSize[1], dataSize[0], dataSize[2]);
    

def inversePermutation(per):
    """Returns the inverse permuation of the permutation per taking axis inversions into account"""
    n = len(per);
    iper = list(per);
    
    #permutation is defined as permuting the axes and then axis inversion
    for i in range(n):
        if per[i] < 0:
            iper[int(abs(per[i])-1)] = -(i + 1);
        else:
            iper[int(abs(per[i])-1)] = (i + 1);
    
    return tuple(iper)

   

def resampleXY(source, sink = None, resolutionSource = (4.0625, 4.0625), resolutionSink = (25, 25), out = sys.stdout):
    """Resample a 2d image"""   
    
    #out.write("Input: %s Output: " % (inputFile, soutputFile))
    data = io.readData(source);
    dataSize = data.shape;
    #print sagittalImageSize;
    
    outputSize = tuple([int(math.ceil(dataSize[i] *  resolutionSource[i]/resolutionSink[i])) for i in range(2)]);
    out.write(("resampleData: Imagesize: %d, %d " % dataSize) + ("Resampled Imagesize: %d, %d" % outputSize))
    
    # note: cv2.resize reverses x-Y axes
    outputData = cv2.resize(data,  (outputSize[1], outputSize[0]));
    #sagittalImage = scipy.misc.imresize(sagittalImage, outputImageSize, interp = 'bilinear'); #normalizes images -> not usefull for stacks !
    out.write("resampleData: resized Image size: %d, %d " % outputData.shape)
    
    return io.writeData(sink, outputData);


def resampleXYParallel(arg):
    """Resampling function to use for parallel resmapling of mage slices"""
    
    inputFile = arg[0];
    outputFile = arg[1];
    resolutionSource = arg[2];
    resolutionSink = arg[3];
    ii = arg[4];
    nn = arg[5];
    
    pw = ProcessWriter(ii);
    pw.write("resampleData: resampling in XY: image %d / %d" % (ii, nn))
    
    data = numpy.squeeze(io.readData(inputFile, z = ii));
    self.resampleXY(data, sink = outputFile, resolutionSource = resolutionSource,  resolutionSink = resolutionSink, out = pw);


def resampleData(source, sink = None, resolutionSource = (4.0625, 4.0625, 3), orientation = None, resolutionSink = (25, 25, 25), processingDirectory = None, processes = 1, cleanup = True):
    """Resample data in source in resolution and orientation"""   
    """Notes: resolutions are assumed to be given for the axes of the intrinsic orientation of the data and reference as when viewed by matplotlib or ImageJ,
              in particular images are loaded with shape (y,x) into numpy arrays
              the scaling is determined by taking into acccount potential reorientation of the data
              orientation: permuation of 1,2,3 with potential sign, indicating which axes map onto the reference axes, 
              negative sign indicates reversal of that particular axes        
    """
    
    #fix (y,x,z) image array representation
    resolutionSource, resolutionSink = self.fixResolutions(resolutionSource, resolutionSink);
    orientation = fixOrientation(orientation);
    
    #orient actual resolutions onto reference resolution
    if not orientation is None:
        resolutionSource = tuple(resolutionSource[int(abs(i))-1] for i in orientation);
    
    if processingDirectory == None:
        processingDirectory = tempfile.mkdtemp();
    
    nZ = io.dataZSize(source);
    
    #resacle in x y
    pool = multiprocessing.Pool(processes=processes);
    argdata = [];
    for i in range(nZ):
        argdata.append( (source, os.path.join(processingDirectory, 'resample_%d.tif' % i), resolutionSource, resolutionSink, i, nZ) );  
    pool.map(self.resampleXYParallel, argdata);
    
    #rescale in z we read 2d images as y,x
    fn = os.path.join(processingDirectory, 'resample_%d.tif' % 0);
    data = io.readData(fn);
    dataSize = data.shape;
    zImage = numpy.zeros((nZ, dataSize[0], dataSize[1]));
    zImage[0,:,:] = data;
    
    for i in range(1,nZ):
        if i % 10 == 0:
            print "resampleData; reading %d/%d" % (i, nZ);
        fn = os.path.join(processingDirectory, 'resample_%d.tif' % i);
        zImage[i,:,:] = io.readData(fn);

    resizedZSize = int(math.ceil(nZ * resolutionSource[2]/resolutionSink[2]));
    resampledData = numpy.zeros((resizedZSize, dataSize[0], dataSize[1]));

    for i in range(dataSize[0]):
        if i % 10 == 0:
            print "resampleData: processing %d/%d" % (i, dataSize[0])
        #resampledImage[:, iImage ,:] =  scipy.misc.imresize(zImage[:,iImage,:], [resizedZAxisSize, sagittalImageSize[1]] , interp = 'bilinear'); 
        #cv2.resize takes reverse order of sizes !
        resampledData[:, i ,:] =  cv2.resize(zImage[:,i,:], (dataSize[1], resizedZSize));
    

    #account for using (z,y,x) array representation -> (y,x,z)
    resampledData = resampledData.transpose([1,2,0]);
    
    print "resampleData: resampled image size:" + str(resampledData.shape)  
    
    if cleanup:
        shutil.rmtree(processingDirectory);

    if not orientation is None:
        
        #reorient
        per = [abs(x)-1 for x in orientation];
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
    
    if sink == []:
        if io.isFileExpression(source):
            sink = os.path.split(source);
            sink = sink[0];
        elif isinstance(source, basestring):
            sink = source + '_resample.tif';
        else:
            raise RuntimeError('resampleData: automatic sink naming not supported for no n str source!');
    
    return io.writeData(sink, resampledData);
    

def resamplePoints(source, dataSize, sink = None, resolutionSource = (4.0625, 4.0625, 3), resolutionSink = (25, 25, 25), orientation = None):
    """Resample Points to map from original data to the coordinates of the resampled image"""
    
    #fix (y,x,z) image array representation
    resolutionSource, resolutionSink = self.fixResolutions(resolutionSource, resolutionSink);
    orientation = self.fixOrientation(orientation);
    
    #orient actual resolutions onto reference resolution
    if not orientation is None:
        resolutionSource = tuple(resolutionSource[int(abs(i))-1] for i in orientation);
    
    points = io.readPoints(source);
    
    #datasize of file pattern
    if isinstance(dataSize, basestring):
        dataSize = io.dataSize(dataSize);
    else:
        dataSize = self.fixDataSize(dataSize);

    #orient resolutions to reference resolutions
    if not orientation is None:
        resolutionSource = tuple(resolutionSource[int(abs(i))-1] for i in orientation);   
    
    #scaling factors
    scale = [int(math.ceil(dataSize[i] *  resolutionSource[i]/resolutionSink[i])) / float(dataSize[i]) for i in range(3)];
    repoints = points.copy();
    for i in range(3):    
        repoints[:,i] = repoints[:,i] * scale[i];
               
    #permute for non trivial orientation
    if not orientation is None:
        per = [int(abs(x))-1 for x in orientation];
        repoints = repoints[:,per];
        
        for i in range(3):
            if orientation[i] < 0:
                repoints[:,i] = dataSize[i] * scale[i] - repoints[:,i];
      
    return io.writePoints(sink, repoints);


     
def resamplePointsInverse(source, dataSize, sink = None, resolutionSource = (4.0625, 4.0625, 3), resolutionSink = (25, 25, 25), orientation = None):
    """Resample Points to map from the coordinates of the resampled image to the original data"""
   
    #fix (y,x,z) image array representation
    resolutionSource, resolutionSink = self.fixResolutions(resolutionSource, resolutionSink);
    orientation = fixOrientation(orientation);
    
    
    #orient actual resolutions onto reference resolution
    if not orientation is None:
        resolutionSource = tuple(resolutionSource[int(abs(i))-1] for i in orientation);    
    
    #datasize of file pattern
    if isinstance(dataSize, basestring):
        dataSize = dataSize(dataSize);
    else:
        dataSize = self.fixDataSize(dataSize);
    
    
    points = io.readPoints(source);
    
    #orient resolutions to reference resolutions
    if not orientation is None:
        resolutionSource = tuple(resolutionSource[int(abs(i))-1] for i in orientation);  

    scale = [int(math.ceil(dataSize[i] *  resolutionSource[i]/resolutionSink[i])) / float(dataSize[i]) for i in range(3)];
    
    rpoints = points.copy();    
    
    #invert axis inversion and permutations    
    if not orientation is None:
        iorientation = self.inversePermutation(orientation);
        #print orientation, iorientation        
        
        #invert permuation
        per = [int(abs(x))-1 for x in iorientation];
        rpoints = rpoints[:,per];
        
        #resampledImage = resampledImage.transpose((0,2,1));
        for i in range(3):
            if iorientation[i] < 0:
                rpoints[:,i] = scale[i] * dataSize[i] - rpoints[:,i];
    
    #scale points
    for i in range(3):   
        rpoints[:,i] = rpoints[:,i] / scale[i];    
    
    return io.writePoints(sink, rpoints);




def test():
    """Test Resampling Module"""
    #import iDISCO.Alignment.Resampling as self
    #reload(self)
    from iDISCO.Parameter import iDISCOPath 
    import iDISCO.IO.IO as io
    import os

    fn = os.path.join(iDISCOPath(), 'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z\d{4}.ome.tif');
    outfn = os.path.join(iDISCOPath(), "Test/Data/Resampling/test.raw")
    
    print "Making subsampled stack " + outfn
    data = self.resampleData(fn, sink = None, resolutionSource = (1,1,1), orientation = (1,2,3), resolutionSink = (10,10,2));
    print data.shape
    io.writeData(outfn, data)   
    

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
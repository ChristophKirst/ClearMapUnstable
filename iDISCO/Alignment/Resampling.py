# -*- coding: utf-8 -*-
"""
Resample data to match Allen Brain Atlas

Notes: imread imsave functions assume (y,x) array representation of image. 
       here all parameters are assumed to be given as x,y,z and are fixed to match the internal array representation
       points are 

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import multiprocessing  
import tempfile
import re
import matplotlib.pyplot as plt
import tifffile as tiff
import os
import math
import numpy
import shutil
import cv2

from iDISCO.Utils.ProcessWriter import ProcessWriter;


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
    
    
def dataSize(imageFilePattern):
    """Determine full size from raw data in (x,y,z) order (not the array format (y,x,z))"""
    
    imageDirectory, listOfImages = self.readFileList(imageFilePattern);
    nz = len(listOfImages);
    
    if nz == 0:
        raise RuntimeError("dataSize: no files match: %s" % imageFilePattern);
    
    imagefile = os.path.join(imageDirectory, listOfImages[0]);
    sagittalImage = plt.imread(imagefile);   
    
    return  (sagittalImage.shape[1], sagittalImage.shape[0], nz)


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


def fixResolutions(resolutionData, resolutionReference):
    """Fix orientation as images are represented as (y,x,z) arrays but the user provides orientation as how to map x,y,z"""

    rd = (resolutionData[1], resolutionData[0], resolutionData[2]);
    rr = (resolutionReference[1], resolutionReference[0], resolutionReference[2]);
    
    return (rd,rr);  
   
def fixDataSize(datasize):
    return (datasize[1], datasize[0], datasize[2]);
    

def inversePermutation(per):
    
    n = len(per);
    iper = list(per);
    
    #permutation is defined as permuting the axes and then axis inversion
    for i in range(n):
        if per[i] < 0:
            iper[int(abs(per[i])-1)] = - (i + 1);
        else:
            iper[int(abs(per[i])-1)] = (i + 1);
    
    return tuple(iper)

   

def resampleXY(inputFile, outputFile, resolutionData = (4.0625, 4.0625), resolutionReference = (25, 25), out = sys.stdout):
    
    #out.write("Input: %s Output: " % (inputFile, soutputFile))
    sagittalImage = plt.imread(inputFile);
    sagittalImageSize = sagittalImage.shape;
    
    outputImageSize = tuple([int(math.ceil(sagittalImageSize[i] *  resolutionData[i]/resolutionReference[i])) for i in range(2)]);
    out.write(("resampleData: Imagesize: %d, %d " % sagittalImageSize) + ("Resampled Imagesize: %d, %d" % outputImageSize))
    
    # note: cv2.resize reverses x-Y axes
    sagittalImage = cv2.resize(sagittalImage,  (outputImageSize[1], outputImageSize[0]));
    #sagittalImage = scipy.misc.imresize(sagittalImage, outputImageSize, interp = 'bilinear'); #normalizes images -> not usefull for stacks !
    out.write("resampleData: resized Image size: %d, %d " % sagittalImage.shape)
    
    tiff.imsave(outputFile, sagittalImage);

def resampleXYParallel(arg):
    """Resampling function to use for parallel resmapling of mage slices"""

    inputFile = arg[0];
    outputFile = arg[1];
    resolutionData = arg[2];
    resolutionReference = arg[3];
    ii = arg[4];
    nn = arg[5];
    
    pw = ProcessWriter(ii);
    pw.write("resampleData: resampling in XY: image %d / %d" % (ii, nn))
    
    self.resampleXY(inputFile, outputFile, resolutionData = resolutionData,  resolutionReference = resolutionReference, out = pw);


def resampleData(imageFilePattern, outputFile = None, resolutionData = (4.0625, 4.0625, 3), orientation = None, resolutionReference = (25, 25, 25), processingDirectory = None, processes = 1, cleanup = True):
    """Resample data to the resolution and orientation of the reference image"""   
    """Notes: resolutions are assumed to be given for the axes of the intrinsic orientation of the data and reference as when viewed by matplotlib or ImageJ,
              in particular images are loaded with shape (y,x) into numpy arrays
              the scaling is determined by taking into acccount potential reorientation of the data
              orientation: permuation of 1,2,3 with potential sign, indicating which axes map onto the reference axes, 
              negative sign indicates reversal of that particular axes        
    """
    
    #fix (y,x,z) image array representation
    resolutionData, resolutionReference = self.fixResolutions(resolutionData, resolutionReference);
    orientation = fixOrientation(orientation);
    
    #orient actual resolutions onto reference resolution
    if not orientation is None:
        resolutionData = tuple(resolutionData[int(abs(i))-1] for i in orientation);
    
    if processingDirectory == None:
        processingDirectory = tempfile.mkdtemp();
    
    imageDirectory, listOfImages = self.readFileList(imageFilePattern);
    nImage = len(listOfImages);
    
    #resacle in x y
    pool = multiprocessing.Pool(processes=processes);
    argdata = [];
    for i in range(nImage):
        argdata.append((os.path.join(imageDirectory, listOfImages[i]), os.path.join(processingDirectory, listOfImages[i] + '.resample.tif'), resolutionData, resolutionReference, i, nImage) );  
    pool.map(self.resampleXYParallel, argdata);
    
    #rescale in z, tiff takes z,y,x coords, we read 2d images as y,x already !
    imageFileName = os.path.join(processingDirectory, listOfImages[0] + '.resample.tif');
    sagittalImage = plt.imread(imageFileName);
    sagittalImageSize = sagittalImage.shape;
    zImage = numpy.zeros((nImage, sagittalImageSize[0], sagittalImageSize[1]));
    for iImage in range(1,nImage):
        if iImage % 10 == 0:
            print "resampleData; reading %d/%d" % (iImage, nImage);
        imageFileName = os.path.join(processingDirectory, listOfImages[iImage] + '.resample.tif');
        zImage[iImage,:,:] = plt.imread(imageFileName);

    resizedZAxisSize = int(math.ceil(nImage * resolutionData[2]/resolutionReference[2]));
    resampledImage = numpy.zeros((resizedZAxisSize, sagittalImageSize[0], sagittalImageSize[1]));

    for iImage in range(sagittalImageSize[0]):
        if iImage % 10 == 0:
            print "resampleData: processing %d/%d" % (iImage, sagittalImageSize[0])
        #resampledImage[:, iImage ,:] =  scipy.misc.imresize(zImage[:,iImage,:], [resizedZAxisSize, sagittalImageSize[1]] , interp = 'bilinear'); 
        #cv2.resize takes reverse order of sizes !
        resampledImage[:, iImage ,:] =  cv2.resize(zImage[:,iImage,:], (sagittalImageSize[1], resizedZAxisSize));
    
    
    print "resampleData: resampled image size:" + str(resampledImage.shape)    
    
    if cleanup:
        shutil.rmtree(processingDirectory);
    
    if not orientation is None:
        #account for using (z,y,x) array representation -> (y,x,z)
        resampledImage = resampledImage.transpose([1,2,0]);

        #reorient
        per = [abs(x)-1 for x in orientation];
        resampledImage = resampledImage.transpose(per);
    
        #reverse orientation after permuting e.g. (-2,1) brings axis 2 to first axis and we can reorder there
        if orientation[0] < 0:
            resampledImage = resampledImage[::-1, :, :];
        if orientation[1] < 0:
            resampledImage = resampledImage[:, ::-1, :]; 
        if orientation[2] < 0:
            resampledImage = resampledImage[:, :, ::-1];
        
        #bring back from y,x,z to z,y,x
        resampledImage = resampledImage.transpose([2,0,1]);
    
    
    if outputFile == None:
        return resampledImage.transpose([1,2,0]); # return as (y,x,z) array 
    else:
        if outputFile == []:
            if imageDirectory[-1] == '/' or imageDirectory[-1] == '\\' :
                outputFile = imageDirectory[0:-1] + '_resample.tif';
            else:
                outputFile = imageDirectory + '_resample.tif';
                
        #print resampledImage.max();
        tiff.imsave(outputFile, resampledImage.astype('uint16'));
        return outputFile;
    

def resamplePoints(points, datasize, resolutionData = (4.0625, 4.0625, 3), resolutionReference = (25, 25, 25), orientation = None):
    """Resample Points to map from original data to the coordinates of the resampled image"""
    
    #fix (y,x,z) image array representation
    resolutionData, resolutionReference = self.fixResolutions(resolutionData, resolutionReference);
    orientation = fixOrientation(orientation);
    
    #orient actual resolutions onto reference resolution
    if not orientation is None:
        resolutionData = tuple(resolutionData[int(abs(i))-1] for i in orientation);
    
    #datasize of file pattern
    if isinstance(datasize, basestring):
        datasize = dataSize(datasize);
    datasize = self.fixDataSize(datasize);

    #orient resolutions to reference resolutions
    if not orientation is None:
        resolutionData = tuple(resolutionData[int(abs(i))-1] for i in orientation);   
    
    #scaling factors
    scale = [int(math.ceil(datasize[i] *  resolutionData[i]/resolutionReference[i])) / float(datasize[i]) for i in range(3)];
    repoints = points.copy();
    for i in range(3):    
        repoints[:,i] = repoints[:,i] * scale[i];
               
    #permute for non trivial orientation
    if not orientation is None:
        per = [int(abs(x))-1 for x in orientation];
        repoints = repoints[:,per];
        
        for i in range(3):
            if orientation[i] < 0:
                repoints[:,i] = datasize[i] * scale[i] - repoints[:,i];
      
    return repoints;


     
def resamplePointsInverse(points, datasize, resolutionData = (4.0625, 4.0625, 3), resolutionReference = (25, 25, 25), orientation = None):
    """Resample Points to map from the coordinates of the resampled image to the original data"""
   
    #fix (y,x,z) image array representation
    resolutionData, resolutionReference = self.fixResolutions(resolutionData, resolutionReference);
    orientation = fixOrientation(orientation);
    
    
    #orient actual resolutions onto reference resolution
    if not orientation is None:
        resolutionData = tuple(resolutionData[int(abs(i))-1] for i in orientation);    
    
    #datasize of file pattern
    if isinstance(datasize, basestring):
        datasize = dataSize(datasize);
    datasize = self.fixDataSize(datasize);
    
    #orient resolutions to reference resolutions
    if not orientation is None:
        resolutionData = tuple(resolutionData[int(abs(i))-1] for i in orientation);  

    scale = [int(math.ceil(datasize[i] *  resolutionData[i]/resolutionReference[i])) / float(datasize[i]) for i in range(3)];
    
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
                rpoints[:,i] = scale[i] * datasize[i] - rpoints[:,i];
    
    #scale points
    for i in range(3):   
        rpoints[:,i] = rpoints[:,i] / scale[i];    
    
    return rpoints;


def test():
    """Resampling"""  
    fn = os.path.split(self.__file__);
    fn = os.path.join(fn[0], '../Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z\d{4}.ome.tif');
    
    print "Searching images of form: " + fn;
    
    fp, fl = self.readFileList(fn);
    print "Found " + str(len(fl)) + " images!"
    
    fp = os.path.split(self.__file__);
    outfn = os.path.join(fp[0], "../Test/Data/Resample/test.tif")
    
    print "Making subsampled stack " + outfn
    
    resampleData(fn, outfn, zResolutionData = 1, zResolutionReference = 2, hemisphere = 1)
        
    

if __name__ == "__main__":
    
    self.test();
    

# -*- coding: utf-8 -*-
"""
Resample data to match Allen Brain Atlas

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];


import multiprocessing  
import tempfile
import scipy.misc
import re
import matplotlib.pyplot as plt
import tifffile as tiff
import os
import math
import numpy
import shutil

import cv2



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

def resampleXY(inputFile, outputFile, xInResolution = 4.0625, yInResolution = 4.0625, xOutResolution = 25, yOutResolution = 25):
    
    #print inputFile
    #print outputFile    
    
    sagittalImage = plt.imread(inputFile);
    sagittalImageSize = sagittalImage.shape;
    #outputImageSize = [int(math.ceil(x)) for x in (sagittalImageSize[0] * xInResolution/xOutResolution, sagittalImageSize[1] * yInResolution/yOutResolution)];
    outputImageSize = (int(math.ceil(sagittalImageSize[1] * yInResolution/yOutResolution)), int(math.ceil(sagittalImageSize[0] * xInResolution/xOutResolution)));
    #print outputImageSize    
    
    #sagittalImage = scipy.misc.imresize(sagittalImage, outputImageSize, interp = 'bilinear');
    sagittalImage = cv2.resize(sagittalImage,  outputImageSize);
    
    tiff.imsave(outputFile, sagittalImage); 


def resmapleXYParallel(arg):
    inputFile = arg[0];
    outputFile = arg[1];
    xInResolution = arg[2];
    yInResolution = arg[3];
    xOutResolution = arg[4];
    yOutResolution = arg[5];
    
    self.resampleXY(inputFile, outputFile, xInResolution, yInResolution, xOutResolution, yOutResolution);


def resampleData(imageFilePattern, outputFile = None, xInResolution = 4.0625, yInResolution = 4.0625, zInResolution = 3, xOutResolution = 25, yOutResolution  = 25, zOutResolution = 25, processingDirectory = None, datasetPrefix = None, processes = 1, cleanup = True, hemisphere = 0):
    if processingDirectory == None:
        processingDirectory = tempfile.mkdtemp();
    
    imageDirectory, listOfImages = self.readFileList(imageFilePattern);
    nImage = len(listOfImages);
    
    #resacle in x y
    pool = multiprocessing.Pool(processes=processes);
    argdata = [];
    for i in range(nImage):
        argdata.append((os.path.join(imageDirectory, listOfImages[i]), os.path.join(processingDirectory, listOfImages[i] + '.resample.tif'), xInResolution, yInResolution, xOutResolution, yOutResolution) );  
    pool.map(self.resmapleXYParallel, argdata);
    
    #rescale in z
    imageFileName = os.path.join(processingDirectory, listOfImages[0] + '.resample.tif');
    sagittalImage = plt.imread(imageFileName);
    sagittalImageSize = sagittalImage.shape;
    zImage = numpy.zeros((nImage, sagittalImageSize[0], sagittalImageSize[1]));
    for iImage in range(1,nImage):
        print "resampleData; reading %d/%d" % (iImage, nImage);
        imageFileName = os.path.join(processingDirectory, listOfImages[iImage] + '.resample.tif');
        zImage[iImage,:,:] = plt.imread(imageFileName);

    resizedZAxisSize = int(math.ceil(nImage * zInResolution/zOutResolution));
    resampledImage = numpy.zeros((resizedZAxisSize, sagittalImageSize[0], sagittalImageSize[1]));

    for iImage in range(sagittalImageSize[0]):
        print "resampleData: processing %d/%d" % (iImage, sagittalImageSize[0])
        resampledImage[:, iImage ,:] =  scipy.misc.imresize(zImage[:,iImage,:], [resizedZAxisSize, sagittalImageSize[1]] , interp = 'bilinear'); 
        resampledImage[:, iImage ,:] =  cv2.resize(zImage[:,iImage,:], (sagittalImageSize[1], resizedZAxisSize,));
    
    if cleanup:
        shutil.rmtree(processingDirectory);
    
    if hemisphere > 0:
        #resampledImage = resampledImage.transpose((0,2,1));
        resampledImage = resampledImage[:,:,::-1];
    
    if outputFile == None:
        return resampledImage;
    else:
        if outputFile == []:
            if imageDirectory[-1] == '/' or imageDirectory[-1] == '\\' :
                outputFile = imageDirectory[0:-1] + '_resample.tif';
            else:
                outputFile = imageDirectory + '_resample.tif';
                
        print resampledImage.max();
        tiff.imsave(outputFile, resampledImage.astype('uint16'));


def resamplePoints(points, datasize, xInResolution = 4.0625, yInResolution = 4.0625, zInResolution = 3, xOutResolution = 25, yOutResolution  = 25, zOutResolution = 25, hemisphere = 0):
    """Resample Points to map fromoriginal data to the coordinates of the resampled image"""
    
    #todo figure out possilble inversion and exchanges of axes
    
    #scaling factors
    
    scale = [int(math.ceil(x) for x in [datasize[0] *  xInResolution/xOutResolution,  datasize[1] * yInResolution/yOutResolution, datasize[2] * zInResolution/zOutResolution]];
    
    
    
def resamplePointsInverse(points,  xInResolution = 4.0625, yInResolution = 4.0625, zInResolution = 3, xOutResolution = 25, yOutResolution  = 25, zOutResolution = 25, hemisphere = 0):
    """Resample Points to map from the coordinates of the resampled image to the original data"""
    




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
    
    resampleData(fn, outfn, zInResolution = 1, zOutResolution = 2, hemisphere = 1)
        
    

if __name__ == "__main__":
    
    self.test();
    

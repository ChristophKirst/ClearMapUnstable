# -*- coding: utf-8 -*-

"""
Spot Detection Algorithm

read data and write points


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];
#self = sys.modules['iDISCO.ImageProcessing.SpotDetection']
#self = sys.modules['__main__'];

import numpy

import cv2
#import mahotas as morph
#from scipy.ndimage.morphology import binary_opening, grey_opening
from scipy.ndimage.filters import correlate
#from scipy.signal import fftconvolve
#from skimage.morphology import reconstruction

#from scipy.ndimage import maximum_filter
from scipy.ndimage.measurements import label, center_of_mass

#import scipy
#from skimage.filter.rank import tophat
#from skimage.measure import regionprops
from mahotas import regmin


#own routines save memory and use faster data types
#from iDISCO.ImageProcessing.Convolution import convolve
from iDISCO.ImageProcessing.GreyReconstruction import reconstruction

from iDISCO.Parameter import ImageProcessingParameter
from iDISCO.ImageProcessing.Filter.StructureElement import structureElement
from iDISCO.ImageProcessing.Filter.FilterKernel import filterKernel

from iDISCO.Utils.Timer import Timer
from iDISCO.Visualization.Plot import plotTiling, plotOverlayLabel

   
def hMaxTransform(img, h):
    """Calculates h maximum transformation of an image."""
    #seed = img.copy();
    #seed[seed < h] = h; # catch errors for uint subtraction !
    #img = img.astype('float16'); # float32 ?
    return reconstruction(img - h, img);


def regionalMax(img):
    """Calculates regional maxima of an image."""
    return regmin(-img);
    
    
def extendedMax(img, h):
    """Calculates extened h maxima of an image."""
    #h max transformimport scipy
    img = hMaxTransform(img, h);
        
    #regional max
    return regionalMax(img);
    
def removeBackground(img, verbose = False, out = sys.stdout, parameter = ImageProcessingParameter()):
    """Remove Background step of Spot Detection Algorithm"""
    timer = Timer();
    
    # background subtraction in each slice
    timer.reset();
    se = structureElement('Disk', parameter.Parameter.Background).astype('uint8');
    for z in range(img.shape[2]):
         #img[:,:,z] = img[:,:,z] - grey_opening(img[:,:,z], structure = structureElement('Disk', (30,30)));
         #img[:,:,z] = img[:,:,z] - morph.grey_opening(img[:,:,z], structure = self.structureELement('Disk', (150,150)));
         img[:,:,z] = img[:,:,z] - cv2.morphologyEx(img[:,:,z], cv2.MORPH_OPEN, se)    
    out.write(timer.elapsedTime(head = 'Background') + '\n');
    
    if verbose:
        plotTiling(10*img);
        
    return img

def dogFilter(img, verbose = False, out = sys.stdout, parameter = ImageProcessingParameter()):
    """Difference of Gaussians filter step of Spot Detection Algorithm"""
    
    timer = Timer();
    
    #DoG filter
    timer.reset();
    fdog = filterKernel(ftype = 'DoG', size = parameter.Parameter.Dog);
    fdog = fdog.astype('float32');
    img = img.astype('float32');
    #img = correlate(img, fdog);
    #img = scipy.signal.correlate(img, fdog);
    img = correlate(img, fdog);
    #img = convolve(img, fdog, mode = 'same');
    img[img < 0] = 0;
    out.write(timer.elapsedTime(head = 'DoG') + '\n');    
    
    if verbose:
        imgp = img;
        imgp[imgp > 200] = 200;
        plotTiling(imgp);
    
    return img

def findExtendedMaxima(img, verbose = False, out = sys.stdout, parameter = ImageProcessingParameter()):
    """Find extended maxima of step in Spot Detection Algorithm"""  
        
    timer = Timer();
    
    # extended maxima
    timer.reset(); 
    imgmax = hMaxTransform(img, parameter.Parameter.HMax);
    imgmax = regionalMax(imgmax);
    imgmax = imgmax.astype('float32') * img;
    th = parameter.Parameter.Threshold;
    imgmax = imgmax > th;
    out.write(timer.elapsedTime(head = 'Extended Max') + '\n');
    
    if verbose:
        #plotTiling(img)
        plotOverlayLabel(img * 0.01, imgmax.astype('int64'), alpha = False);
        #plotOverlayLabel(img, imgmax.astype('int64'), alpha = True)     
        
    return imgmax


def findCenterOfMaxima(img, imgmax, verbose = False, out = sys.stdout):
    """Find center of the maxima step in Spot Detection Algorithm"""  
    
    timer = Timer();    
    
    #center of maxima
    timer.reset();
    imglab, nlab = label(imgmax);
    
    if nlab > 0:
        centers = numpy.array(center_of_mass(img, imglab, index = numpy.arange(1, nlab)));    
        out.write(timer.elapsedTime(head = 'Cell Centers') + '\n');
    
        if verbose:        
            #plotOverlayLabel(img * 0.01, imglab, alpha = False);
            #plotTiling(img)
            imgc = numpy.zeros(img.shape);
            for i in range(centers.shape[0]):
                imgc[centers[i,0], centers[i,1], centers[i,2]] = 1;
            plotOverlayLabel(img, imgc, alpha = False);
            #plotOverlayLabel(img, imgmax.astype('int64'), alpha = True)     
    
        #return centers, imglab, mask
        cintensity = numpy.array([img[centers[i,0], centers[i,1], centers[i,2]] for i in range(centers.shape[0])]);        
    
        return ( centers, cintensity );
    else:
        out.write(timer.elapsedTime(head = 'Cell Centers'));
        out.wrtie('No Cells found !');
        return ( numpy.zeros((0,3)), numpy.zeros(0) );



def detectCells(img, verbose = False, out = sys.stdout, parameter = ImageProcessingParameter()):
    """Detect Cells import scipyin 3d grayscale img using DoG filtering and maxima dtection"""
    # effectively removeBackground, dogFilter, 
    # we do processing steps in place to save memory    
    
    timer = Timer();
    
    # normalize data -> to check
    #img = img.astype('float');
    #dmax = 0.075 * 65535;
    #ids = img > dmax;
    #img[ids] = dmax;
    #img /= dmax; 
    #out.write(timer.elapsedTime(head = 'Normalization'));
    #img = dataset[600:1000,1600:1800,800:830];
    #img = dataset[600:1000,:,800:830];
    
    
    # background subtraction in each slice
    timer.reset();
    se = structureElement('Disk', parameter.Parameter.Background).astype('uint8');
    for z in range(img.shape[2]):
         #img[:,:,z] = img[:,:,z] - grey_opening(img[:,:,z], structure = structureElement('Disk', (30,30)));
         #img[:,:,z] = img[:,:,z] - morph.grey_opening(img[:,:,z], structure = self.structureELement('Disk', (150,150)));
         img[:,:,z] = img[:,:,z] - cv2.morphologyEx(img[:,:,z], cv2.MORPH_OPEN, se)    
    out.write(timer.elapsedTime(head = 'Background'));
    
    if verbose:
        plotTiling(10*img, maxtiles = 9);    
    
    # mask
    #timer.reset();
    #if mask == None: #explicit mask
    #    mask = img > 0.01;
    #    mask = binary_opening(mask, self.structureELement('Disk', (3,3,3)));
    #img[img < 0.01] = 0; # masking in place  # extended maxima
    
    
    #out.write(timer.elapsedTime(head = 'Mask'));
    
    
    #DoG filter
    timer.reset();
    fdog = filterKernel(ftype = 'DoG', size = parameter.Parameter.Dog);
    fdog = fdog.astype('float32');
    img = img.astype('float32');
    #img = correlate(img, fdog);
    #img = scipy.signal.correlate(img, fdog);
    img = correlate(img, fdog);
    #img = convolve(img, fdog, mode = 'same');
    img[img < 0] = 0;
    out.write(timer.elapsedTime(head = 'DoG'));
    
    if verbose:
        imgp = img;
        imgp[imgp > 200] = 200;
        plotTiling(imgp[:,:,10:20], maxtiles = 9);
    
#    imax = img.max();
#    if imax == 0:
#        imax = 1;
#    img /= imax;
    
    # extended maxima
    timer.reset(); 
    imgmax = hMaxTransform(img, parameter.Parameter.HMax);
    imgmax = regionalMax(imgmax);
    imgmax = imgmax.astype('float32') * img;
    th = parameter.Parameter.Threshold;
    imgmax = imgmax > th;
    out.write(timer.elapsedTime(head = 'Extended Max'));
    
    if verbose:
        #plotTiling(img)
        plotOverlayLabel(img * 0.01, imgmax.astype('int64'), alpha = False);
        #plotOverlayLabel(img, imgmax.astype('int64'), alpha = True)     
    
    
    #center of maxima
    timer.reset();
    imglab, nlab = label(imgmax);
    
    if nlab > 0:
        centers = numpy.array(center_of_mass(img, imglab, index = numpy.arange(1, nlab)));    
        out.write(timer.elapsedTime(head = 'Cell Centers'));
    
        if verbose:        
            plotOverlayLabel(img * 0.01, imglab, alpha = False);
            #plotTiling(img)
            imgc = numpy.zeros(img.shape);
            for i in range(centers.shape[0]):
                imgc[centers[i,0], centers[i,1], centers[i,2]] = 1;
            plotOverlayLabel(img * 0.01, imgc, alpha = False);
            #plotOverlayLabel(img, imgmax.astype('int64'), alpha = True)     
    
        #return centers, imglab, mask
        cintensity = numpy.array([img[centers[i,0], centers[i,1], centers[i,2]] for i in range(centers.shape[0])]);        
    
        return ( centers, cintensity );
    else:
        out.write(timer.elapsedTime(head = 'Cell Centers'));
        out.write('No Cells found !');
        return ( numpy.zeros((0,3)), numpy.zeros(0) );







def test():
    import iDISCO.IO.Imaris as io  
    
    fn = '/home/ckirst/Science/Projects/BrainActivityMap/Data/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims';
    #fn = '/run/media/ckirst/ChristophsBackuk4TB/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims';
    #fn = '/home/nicolas/Windows/Nico/cfosRegistrations/Adult cfos C row 20HF 150524 - Copy.ims';
    #fn = '/home/ckirst/Science/Projects/BrainActivityMap/iDISCO_2015_04/test for spots added spot.ims'
    
    
    
    f = io.openFile(fn);
    dataset = io.readData(f, resolution=0);
    #img = dataset[0:500,0:500,1000:1008];
    #img = dataset[600:1000,1600:1800,800:830];
    img = dataset[500:1500,500:1500,800:809];    
    f.close();
    
    print self
    
    #m = sys.modules['iDISCO.ImageProcessing.SpotDetection']
    c = detectCells(img);
    
    print 'done, found %d cells !' % c.shape[0]


if __name__ == '__main__':
    test();
    
    
    
    

# -*- coding: utf-8 -*-

"""
Maxima Detection Algorithms

Notes:
    -used for finding maxima and cell intensities

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import numpy


#from scipy.ndimage import maximum_filter
from scipy.ndimage.measurements import label, center_of_mass

#import scipy
#from skimage.filter.rank import tophat
#from skimage.measure import regionprops
from mahotas import regmin


#own routines save memory and use faster data types
#from iDISCO.ImageProcessing.Convolution import convolve
from iDISCO.ImageProcessing.GreyReconstruction import reconstruction
from iDISCO.ImageProcessing.Filter.StructureElement import  structureElementOffsets


from iDISCO.Utils.Timer import Timer
from iDISCO.Utils.ParameterTools import writeParameter;
from iDISCO.Visualization.Plot import plotOverlayLabel


##############################################################################
# Basic Transforms
##############################################################################

   
def hMaxTransform(img, h):
    """Calculates h maximum transformation of an image."""
    #seed = img.copy();
    #seed[seed < h] = h; # catch errors for uint subtraction !
    #img = img.astype('float16'); # float32 ?
    return reconstruction(img - h, img);


def regionalMax(img):
    """Calculates regional maxima of an image."""
    return regmin(-img);
    
    
def extendedMax(img, hMax = 0):
    """Calculates extened h maxima of an image."""

    #h max transformimport scipy
    img = hMaxTransform(img, hMax);
        
    #regional max
    return regionalMax(img);



def findExtendedMaxima(img, hMax = 10, threshold = 0, verbose = False, out = sys.stdout, **parameter):
    """Find extended maxima of step in Spot Detection Algorithm"""  
    timer = Timer();
    writeParameter(out = out, head = 'Extended Max:', hMax = hMax, threshold = threshold);
    
    # extended maxima
    timer.reset(); 
    imgmax = hMaxTransform(img, hMax);
    imgmax = regionalMax(imgmax);
    imgmax = imgmax.astype('float32') * img;
    imgmax = imgmax > threshold;
    
    if verbose:
        #plotTiling(img)
        plotOverlayLabel(img * 0.01, imgmax.astype('int64'), alpha = False);
        #plotOverlayLabel(img, imgmax.astype('int64'), alpha = True)     
    
    out.write(timer.elapsedTime(head = 'Extended Max') + '\n');
    return imgmax



def findCenterOfMaxima(img, imgmax, verbose = False, out = sys.stdout, **parameter):
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
        #cintensity = numpy.array([img[centers[i,0], centers[i,1], centers[i,2]] for i in range(centers.shape[0])]);
        
        out.write(timer.elapsedTime(head = 'Cell Centers'));
        #return ( centers, cintensity );
        return centers;
        
    else:
        out.wrtie('Cell Centers: No Cells found !');
        #return ( numpy.zeros((0,3)), numpy.zeros(0) );
        return numpy.zeros((0,3));



def findIntensity(img, centers, intensityMethod = 'Max', intensitySize = (3,3,3), verbose = False, out = sys.stdout, **parameter):
    """Find instensity value around centers in the image img using a box"""
    timer = Timer();    
    writeParameter(out = out, head = 'Cell Intensities:', intensityMethod = intensityMethod, intensitySize = intensitySize);
    
    if centers.shape[0] == 0:
        return numpy.zeros(0);
    
    if intensityMethod is None:
            return numpy.array([img[centers[i,0], centers[i,1], centers[i,2]] for i in range(centers.shape[0])]);        
    
    isize = img.shape;
    #print isize
    
    offs = structureElementOffsets(intensitySize);
    
    if isinstance(intensityMethod, basestring):
        method = eval('numpy.' + intensityMethod.lower());


    intensities = numpy.zeros(centers.shape[0], dtype = img.dtype);
    
    for c in range(centers.shape[0]):
        xmin = int(-offs[0,0] + centers[c,0]);
        if xmin < 0:
            xmin = 0;       
        xmax = int(offs[0,1] + centers[c,0]);
        if xmax > isize[0]:
            xmax = isize[0];
            
        ymin = int(-offs[1,0] + centers[c,1]);
        if ymin < 0:
            ymin = 0;       
        ymax = int(offs[1,1] + centers[c,1]);
        if ymax > isize[1]:
            ymax = isize[1];
            
        zmin = int(-offs[2,0] + centers[c,2]);
        if zmin < 0:
            zmin = 0;       
        zmax = int(offs[1,1] + centers[c,2]);
        if zmax > isize[2]:
            zmax = isize[2];
        
        #print xmin, xmax, ymin, ymax, zmin, zmax
        data = img[xmin:xmax, ymin:ymax, zmin:zmax];
        
        intensities[c] = method(data);
    
    out.write(timer.elapsedTime(head = 'Cell Intensities'));
    return intensities;

    

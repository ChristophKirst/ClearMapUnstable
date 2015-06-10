# -*- coding: utf-8 -*-

"""
Spot Detection Algorithm

read data and write points


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];
self = sys.modules['iDISCO.ImageProcessing.SpotDetection']

import numpy

#import mahotas as morph
#from scipy.ndimage.morphology import binary_opening, grey_opening
from scipy.ndimage.filters import correlate
#from scipy.ndimage import maximum_filter
from scipy.ndimage.measurements import label, center_of_mass

from skimage.morphology import reconstruction

import cv2

#import scipy

#from skimage.filter.rank import tophat
#from skimage.measure import regionprops

from mahotas import regmin

from iDISCO.Utils.Timer import Timer
from iDISCO.ImageProcessing.Filter.StructureElement import structureElement
from iDISCO.ImageProcessing.Filter.FilterKernel import filterKernel
from iDISCO.Visualization.Plot import plotTiling, plotOverlayLabel
   
def hMaxTransform(img, h):
    """Calculates h maximum transformation of an image."""
    #seed = img.copy();
    #seed[seed < h] = h; # catch errors for uint subtraction !
    img = img.astype('float');
    return reconstruction(img - h, img);


def regionalMax(img):
    """Calculates regional maxima of an image."""
    return regmin(-img);
    
    
def extendedMax(img, h):
    """Calculates extened h maxima of an image."""
    #h max transformimport scipy
    img = self.hMaxTransform(img, h);
        
    #regional max
    return self.regionalMax(img);
            


def detectCells(img, verbose = False, out = sys.stdout):
    """Detect Cells import scipyin 3d grayscale img using DoG filtering and maxima dtection"""
    
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
    se = structureElement('Disk', (15,15)).astype('uint8');
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
    fdog = filterKernel(ftype = 'DoG', size = [11, 7, 7]);
    img = correlate(img.astype('float'), fdog);
    #img = scipy.signal.correlate(img, fdog);
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
    imgmax = self.hMaxTransform(img, 20);
    imgmax = self.regionalMax(imgmax);
    imgmax = imgmax.astype('float') * img;
    th = 20;
    imgmax = imgmax > th;
    out.write(timer.elapsedTime(head = 'Extened Max'));
    
    if verbose:
        #plotTiling(img)
        plotOverlayLabel(img * 0.01, imgmax.astype('int64'), alpha = False);
        #plotOverlayLabel(img, imgmax.astype('int64'), alpha = True)     
    
    
    #center of maxima
    timer.reset();
    imglab, nlab = label(imgmax);
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
    
    return (centers, cintensity);
    

    
"""plotTiling(numpy.dstack((10 * img, 1 * (img - res2))))



# run segmentation
img = data[:,:,0:30];

import iDISCO.Segmentation.SpotDetection
import sys
self = sys.modules['iDISCO.Segmentation.SpotDetection'];



function hdImage = hdTransform2(I, h, n)
    disp('Enter into h-d transform');
    smallNumber = 1e-7;
    se = strel(ones(n));
    J = I-h;
    flag = 0;
%      [R C] = size(I);
%    count = 0;
    pI = J*0;
    while flag==0
        tempJ = min( imdilate(J, se), I);
        resud = sum( sum(abs(J-tempJ) ));
        J =  tempJ;
        pI = max(pI,J);
%        count = count +1
        if resud < smallNumber
            flag = 1;
        end
    end

    hdImage = I-pI;
%         figure
%     subplot(2,2,1)
%     imshow(I,[])
%     subplot(2,2,2)
%      imshow(pI,[])
%         subplot(2,2,3)
%      imshow( hdImage,[])
%      colormap('jet')
%      linkaxes;
    disp('Exit from h-d transform');
end 



    
    
import numpy as np
from skimage.morphology import reconstruction

x = np.linspace(0, 4 * np.pi)
y_mask = np.cos(x)

y_seed = y_mask.min() * np.ones_like(x)
y_seed[0] = 0.5
y_seed[-1] = 0
y_rec = reconstruction(y_seed, y_mask)
"""

    
    
    

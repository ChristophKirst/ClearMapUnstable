# -*- coding: utf-8 -*-
"""
Illumination Correction

Notes:
    -  corrects image via (img - bkg) / (flt - bkg) 

Reference: 
    - Fundamentals of Light Microscopy and Electronic Imaging, p 421
    
@author: ckirst
"""


import sys
import numpy
import os

import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

import iDISCO.IO.IO as io

from iDISCO.Utils.Timer import Timer
from iDISCO.Utils.ParameterTools import writeParameter

from iDISCO.Visualization.Plot import plotTiling

from iDISCO.Settings import IDISCOPath



def correctIllumination(img, flatfield = True, background = None, illuminationScaling = True, subStack = None, verbose = False, correctedIlluminationFile = None, out = sys.stdout, **parameter):
    """Correct Illumination"""
    timer = Timer();
    if flatfield is None or isinstance(flatfield, str) or flatfield is True:
        fld = flatfield;
    else:
        fld = "image of size %s" % str(flatfield.shape);
    
    if background is None or isinstance(background, str):
        bkg = background;
    else:
        bkg = "image of size %s" % str(background.shape);
    
    writeParameter(out = out, head = 'Illumination:', flatfield = fld,  background = bkg, illuminationScaling = illuminationScaling);
 
    if not subStack is None:
        x = subStack["x"];
        y = subStack["y"];
    else:
        x = all;
        y = all;
    
    #print "sizes", x, y, img.shape
 
    #read data  
    if flatfield is None:
        return img;
        
    elif flatfield is True:
        # default flatfield correction
        flatfield = flatfieldFromLine(os.path.join(IDISCOPath, "Data/lightsheet_flatfield_correction.csv"), img.shape[0]);
            
    elif isinstance(flatfield, str):
        # 
        if io.isPointFile(flatfield):
            flatfield = flatfieldFromLine(flatfield, img.shape[0]);
        else:
            flatfield = io.readData(flatfield);
    
    ffmean = flatfield.mean();    
    ffmax = flatfield.max();

    #correct for subset
    flatfield = io.readData(flatfield, x = x, y = y);   
    
    background = io.readData(background, x = x, y = y);
    
    if flatfield.shape != img[:,:,0].shape:
        raise RuntimeError("correctIllumination: flatfield does not match image size: %s vs %s" % (flatfield.shape,  img[:,:,0].shape));
    
    #convert to float for scaling
    dtype = img.dtype;
    img = img.astype('float32');
    flatfield = flatfield.astype('float32');
    
    # illumination correction in each slice
    if background is None:
        for z in range(img.shape[2]):
            img[:,:,z] = img[:,:,z] / flatfield;
    else:
        if background.shape != flatfield.shape:
            raise RuntimeError("correctIllumination: background does not match image size: %s vs %s" % (background.shape,  img[:,:,0].shape));        
        background = background.astype('float32');

        flatfield = (flatfield - background);
        for z in range(img.shape[2]):
            img[:,:,z] = (img[:,:,z] - background) / flatfield;
    
        
    # rescale
    if illuminationScaling == "mean":
        # scale back by average flat field correction:
        sf = ffmean;
    elif illuminationScaling == "max":
        sf = ffmax;
    else:
        sf = illuminationScaling;
        
    writeParameter(out = out, head = 'Illumination:',  illuminationScaling = sf);
    
    if not sf is None:
        img = img * sf;
        img = img.astype(dtype);
    
    
    #write result for inspection
    if not correctedIlluminationFile is None:
        if not subStack is None:
            ii = subStack["zSubStackCenterIndices"][0];
            ee = subStack["zSubStackCenterIndices"][1];
            si = subStack["zCenterIndices"][0];
        else:
            si = 0;
            ii = 0;
            ee = -1;
        
        io.writeData(correctedIlluminationFile, img[:,:,ii:ee], startIndex = si );     
    
        
    #plot result for inspection
    if verbose:
        plotTiling(img);

    out.write(timer.elapsedTime(head = 'Illumination') + '\n');    
    return img 
    


def flatfieldFromLine(line, xsize):
    """Creates Flatfield image from a 1d line of intensities"""
    line = io.readPoints(line);

    flatfield = numpy.zeros((xsize, line.size));
    for i in range(xsize):
        flatfield[i,:] = line;
    
    return flatfield;



def flatfieldLineFromRegression(data, sink = None, method = 'polynomial', reverse = None, verbose = False):
    """Create Flatfield line fit from a list of positions and intensities"""
    data = io.readPoints(data);

    # split data
    x = data[:,0]
    y = data[:,1:-1];
    ym = numpy.mean(y, axis = 1);

    if verbose:
        plt.figure()
        for i in range(1,data.shape[1]):
            plt.plot(x, data[:,i]);
        plt.plot(x, ym, 'k');
    
    
    if method == 'polynomial':
        ## fit r^6
        mean = sum(ym * x)/sum(ym)

        def f(x,m,a,b,c,d):
            return a + b * (x-m)**2 + c * (x-m)**4 + d * (x-m)**6;
        
        popt, pcov = curve_fit(f, x, ym, p0 = (mean, 1, 1, 1, .1));
        m = popt[0]; a = popt[1]; b = popt[2];
        c = popt[3]; d = popt[4];

        if verbose:
            print "polynomial fit: %f + %f (x- %f)^2 + %f (x- %f)^4 + %f (x- %f)^6" % (a, b, m, c, m, d, m);

        def fopt(x):
            return f(x, m = m, a = a, b = b, c = c, d = d);
        
        flt = map(fopt, range(0, int(x[-1])));
    
    else: 
        ## Gaussian fit
        
        mean = sum(ym * x)/sum(ym)
        sigma = sum(ym * (x-mean)**2)/(sum(ym))
        
        
        def f(x, a, m, s, b):
            return a * numpy.exp(- (x - m)**2 / 2 / s) + b;
            
        
        popt, pcov = curve_fit(f, x, ym, p0 = (1000, mean, sigma, 400));
        a = popt[0]; m = popt[1]; s = popt[2]; b = popt[3];

        if verbose:
            print "Gaussian fit: %f exp(- (x- %f)^2 / (2 %f)) + %f" % (a, m, s, b);

        def fopt(x):
            return f(x, a = a, m = m, s = s, b = b);
    
    
    if verbose:
        plt.plot(x, flt);
    
    
    if reverse:
        flt.reverse();
    
    
    if sink is None:
        return flt;
    else:
        return io.writePoints(sink, flt);
        



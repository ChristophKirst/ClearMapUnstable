# -*- coding: utf-8 -*-

"""
Spot Detection Algorithm

read data and write points


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

from scipy.ndimage.filters import correlate
#from scipy.signal import fftconvolve

from iDISCO.ImageProcessing.Filter.FilterKernel import filterKernel
from iDISCO.ImageProcessing.BackgroundRemoval import removeBackground
from iDISCO.ImageProcessing.MaximaDetection import findExtendedMaxima, findCenterOfMaxima, findIntensity

from iDISCO.Utils.Timer import Timer
from iDISCO.Utils.ParameterTools import writeParameter
from iDISCO.Visualization.Plot import plotTiling


##############################################################################
# Image processing steps in spot detection
##############################################################################


def dogFilter(img, dogSize = (7,7,11), verbose = False, out = sys.stdout, **parameter):
    """Difference of Gaussians filter step of Spot Detection Algorithm"""
    timer = Timer();
    writeParameter(out = out, head = 'DoG:', dogSize = dogSize);
    
    #DoG filter
    timer.reset();
    fdog = filterKernel(ftype = 'DoG', size = dogSize);
    fdog = fdog.astype('float32');
    img = img.astype('float32');
    #img = correlate(img, fdog);
    #img = scipy.signal.correlate(img, fdog);
    img = correlate(img, fdog);
    #img = convolve(img, fdog, mode = 'same');
    img[img < 0] = 0;
    
    if verbose:
        plotTiling(img);
        
    out.write(timer.elapsedTime(head = 'DoG') + '\n');    
    return img


##############################################################################
# Spot detection
##############################################################################

def detectCells(img, verbose = False, out = sys.stdout, **parameter):
    """Detect Cells in 3d grayscale image using DoG filtering and maxima detection"""
    # effectively removeBackground ->  dogFilter -> hHmax -> trheshold
    # processing steps are done in place to save memory  
    
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
    img = removeBackground(img, verbose = verbose, out = out, **parameter)   
    
    # mask
    #timer.reset();
    #if mask == None: #explicit mask
    #    mask = img > 0.01;
    #    mask = binary_opening(mask, self.structureELement('Disk', (3,3,3)));
    #img[img < 0.01] = 0; # masking in place  # extended maxima
    #out.write(timer.elapsedTime(head = 'Mask'));    
    
    #DoG filter
    img = dogFilter(img, verbose = verbose, out = out, **parameter);
    
    # normalize    
    #    imax = img.max();
    #    if imax == 0:
    #        imax = 1;
    #    img /= imax;
    
    # extended maxima
    imgmax = findExtendedMaxima(img, verbose = verbose, out = out, **parameter);
    
    #center of maxima
    centers = findCenterOfMaxima(img, imgmax, verbose = verbose, out = out);
    
    #intensity of cells
    cintensity = findIntensity(img, centers, verbose = verbose, out = out, **parameter);
    
    return ( centers, cintensity );
        




def test():
    """Test Spot Detection Module"""
    import iDISCO.ImageProcessing.SpotDetection as self
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
    

    #m = sys.modules['iDISCO.ImageProcessing.SpotDetection']
    c = detectCells(img);
    
    print 'done, found %d cells !' % c.shape[0]


    #test intensities:
    import numpy;
    x = numpy.random.rand(30,30,10);
    centers = numpy.array([[0,0,0], [29,29,9]]);
    i = self.findIntensity(x, centers, boxSize = (1,1,1));
    print i

if __name__ == '__main__':
    test();
    
    
    
    

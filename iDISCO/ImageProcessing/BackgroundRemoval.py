# -*- coding: utf-8 -*-
"""
Background Removal

Notes:
    - use morphological opening for back ground removal
    
Created on Fri Aug 14 16:48:27 2015

@author: ckirst
"""

import sys

import cv2 

from iDISCO.ImageProcessing.Filter.StructureElement import structureElement

from iDISCO.Utils.Timer import Timer
from iDISCO.Utils.ParameterTools import writeParameter

from iDISCO.Visualization.Plot import plotTiling


def removeBackground(img, backgroundSize = (15,15), verbose = False, out = sys.stdout, **parameter):
    """Remove Background step of Spot Detection Algorithm"""
    timer = Timer();
    writeParameter(out = out, head = 'Background:', backgroundSize = backgroundSize);
    
    # background subtraction in each slice
    se = structureElement('Disk', backgroundSize).astype('uint8');
    for z in range(img.shape[2]):
         #img[:,:,z] = img[:,:,z] - grey_opening(img[:,:,z], structure = structureElement('Disk', (30,30)));
         #img[:,:,z] = img[:,:,z] - morph.grey_opening(img[:,:,z], structure = self.structureELement('Disk', (150,150)));
         img[:,:,z] = img[:,:,z] - cv2.morphologyEx(img[:,:,z], cv2.MORPH_OPEN, se)
        
    if verbose:
        plotTiling(10*img);

    out.write(timer.elapsedTime(head = 'Background') + '\n');    
    return img 

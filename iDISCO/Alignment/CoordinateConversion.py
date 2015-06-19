# -*- coding: utf-8 -*-
"""
Convert Coordinates to various formats

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

from iDISCO.Alignment.ElastixParameter import ElastixParameter

def scaleCoordinates(points, scale):
    """Scale coordinates"""
    if len(scale) == 1:
        scale = (scale, scale, scale);
    
    for i in range(3):
        points[:,i] = points[:,i] * scale[i];

def transformToImaris(points, scale = (4.0625, 4.0625, 3)):
    """Transform pixel coordinates of cell centers to work in Imaris"""
    return self.scaleCoordinates(points, scale = scale);

def transfromByEalstix(points, elastixParameter = ElastixParameter()):#
    """Transform pixel coordinates of cell centers using a calculated transformation obtained fomr elastix"""
        
    print 'todo'
    






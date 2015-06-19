# -*- coding: utf-8 -*-
"""
Simple Interface to writing VTK files

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

from evtk.hl import pointsToVTK 

def writePoints(filename, points, data = None):
    x = points[:,0];
    y = points[:,1];
    z = points[:,2];
    
    if data == None:
        pointsToVTK(filename, x, y, z);
    else:
        pointsToVTK(filename, x, y, z, data = data)
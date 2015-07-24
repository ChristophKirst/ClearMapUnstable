# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 14:50:57 2015

@author: ckirst
"""

#cimport cython 

import numpy as np
cimport numpy as np

#IM_TYPE = np.int_
#ctypedef np.int_t IM_TYPE_t
#ctypedef np.float_t PT_TYPE_t

#@cython.boundscheck(False)  # turn of bounds-checking for entire function
def voxelizeSphere(np.ndarray[np.float_t, ndim=2] points, int xsize, int ysize, int zsize, float xdiam, float ydiam, float zdiam):
    
    cdef np.ndarray[np.int_t, ndim =3] voximg = np.zeros([xsize, ysize, zsize], dtype=np.int)
    cdef int x, y, z

    cdef int iCentroid = 0;
    cdef int nCentroid = points.shape[0];
    cdef int nSphereIndices = int(xdiam * ydiam * zdiam);

    # precompute indices centered at 0,0,0
    cdef np.ndarray[np.int_t, ndim = 1] xs = np.zeros([nSphereIndices], dtype=np.int);
    cdef np.ndarray[np.int_t, ndim = 1] ys = np.zeros([nSphereIndices], dtype=np.int);
    cdef np.ndarray[np.int_t, ndim = 1] zs = np.zeros([nSphereIndices], dtype=np.int);
    cdef int ns = 0;

    cdef float xdiam2 = xdiam * xdiam / 4;
    cdef float ydiam2 = ydiam * ydiam / 4;
    cdef float zdiam2 = zdiam * zdiam / 4;
    
    for x in range(int(-xdiam/2 + 1), int(xdiam/2 + 1)):
        for y in range(int(-ydiam/2 + 1), int(ydiam/2 + 1)):
            for z in range(int(-zdiam/2 + 1), int(zdiam/2 + 1)):
                if x*x / xdiam2 + y*y / ydiam2 + z*z / zdiam2 < 1:
                    xs[ns] = x; ys[ns] = y; zs[ns] = z;
                    ns += 1;
                    
    cdef int iss = 0;
    cdef float cx0;
    cdef float cy0;
    cdef float cz0;
    
    cdef int cx; 
    cdef int cy;
    cdef int cz;
                    
    for iCentroid in range(nCentroid):
        if ((iCentroid % 25000) == 0):
            print "\nProcessed %d/%d\n" % (iCentroid, nCentroid);
    
        cx0 = points[iCentroid, 0];
        cy0 = points[iCentroid, 1];
        cz0 = points[iCentroid, 2];
        
        for iss in range(ns):
            cx = int(cx0 + xs[iss]);
            cy = int(cy0 + ys[iss]);
            cz = int(cz0 + zs[iss]);
            
            if cx >= 0 and cx < xsize:
                if cy >= 0 and cy < ysize:
                    if cz >= 0 and cz < zsize:
                        voximg[cx,cy,cz] = voximg[cx,cy,cz] + 1;
                        
    return voximg;
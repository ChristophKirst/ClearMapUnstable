# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 13:55:06 2015

@author: kannanuv
"""

#import iDISCO
import numpy
import math
import random
from PIL import Image

# if the small loop bounds are out of image bounds then set them within image bounds
def sanitizeLoopMinAndMax (iLoopMin, iLoopMax, iImageMin, iImageMax):
    if (iLoopMin < iImageMin):
        iLoopMin = iImageMin;    
        
    if (iLoopMax < iImageMin):
        iLoopMax = iImageMin;
        
    if (iLoopMin > iImageMax):
        iLoopMin = iImageMax;
        
    if (iLoopMax > iImageMax):
        iLoopMax = iImageMax;
        
    return iLoopMin, iLoopMax

# write voxelzied image
def writeVoxelizedImage (voxelizedImageFileName, pointList, parameter):
    
    #swap x and y axis
    xVec = pointList[:,1];
    yVec = pointList[:,0];
    zVec = pointList[:,2];
    
    distance = 2;

    #xScale = 2.5; Use this incase the cell detected are still in image pixel space and not real space
    #yScale = 2.5;
    #zScale = 3.0;

    xScale = 1.0;
    yScale = 1.0;
    zScale = 1.0;
    
    nCell = xVec.size;
            
    #Image dimensions
    nX = 450;
    nY = 650;
    nZ = 300;
    
    voxelNature = 0; # 0 for spherical overlapping and 1 for non-overlapping rectangular voxels
    
    xImageMin = 0;
    xImageMax = 450;
    yImageMin = 0;
    yImageMax = 650;
    zImageMin = 0;
    zImageMax = 300;
    
    #Spacing
    xSpacing = (xImageMax - xImageMin) / nX;
    ySpacing = (yImageMax - yImageMin) / nY;
    zSpacing = (zImageMax - zImageMin) / nZ;
    
    voxelImage = numpy.zeros([nX,nY, nZ], float)

    for iCell in range(0, nCell):
        #Update the user of processing every 25000 cells
        if ((iCell % 20) == 0):
            print (str(iCell) + " / " + str(nCell) + " processing completed");
            
        x = xVec[iCell];
        y = yVec[iCell];
        z = zVec[iCell];
        
        #find the smallest volume to be looped
        xLoopMin = int (math.floor ((xVec[iCell] - distance/xScale - xImageMin) / xSpacing));
        yLoopMin = int (math.floor ((yVec[iCell] - distance/yScale - yImageMin) / ySpacing));
        zLoopMin = int (math.floor ((zVec[iCell] - distance/zScale - zImageMin) / zSpacing));
        xLoopMax = int (math.ceil ((xVec[iCell] + distance/xScale - xImageMin) / xSpacing));
        yLoopMax = int (math.ceil ((yVec[iCell] + distance/yScale - yImageMin) / ySpacing));
        zLoopMax = int (math.ceil ((zVec[iCell] + distance/zScale - zImageMin) / zSpacing));
        
        xLoopMin, xLoopMax = sanitizeLoopMinAndMax (xLoopMin, xLoopMax, xImageMin, xImageMax);
        yLoopMin, yLoopMax = sanitizeLoopMinAndMax (yLoopMin, yLoopMax, yImageMin, yImageMax);
        zLoopMin, zLoopMax = sanitizeLoopMinAndMax (zLoopMin, zLoopMax, zImageMin, zImageMax);
        
        for iX in range (xLoopMin, xLoopMax):
            curX = xImageMin + iX * xSpacing;
            for iY in range (yLoopMin, yLoopMax):
                curY = yImageMin + iY * ySpacing;
                for iZ in range (zLoopMin, zLoopMax):
                    curZ = zImageMin + iZ * zSpacing;
                    if (voxelNature == 0): #spherical
                        l2x = math.pow(xScale,2) * math.pow((curX - x),2);
                        l2y = math.pow(yScale,2) * math.pow((curY - y),2);
                        l2z = math.pow(zScale,2) * math.pow((curZ - z),2);
                        l2 = math.sqrt(l2x + l2y + l2z);
                        if (l2 < distance):
                            voxelImage[iX, iY, iZ] = voxelImage[iX, iY, iZ] + 1;
                            
    resultImage = Image.fromarray((voxelImage).astype(numpy.uint16))
    resultImage.save('voxelizedImage.tif')
        
# Test program for the voxelization        
nPoint = 100;
pointList = numpy.zeros([nPoint,3], float);
for iPoint in range(0, nPoint):
    for iCoord in range(0, 3):
        if (iCoord == 0):
            pointList[iPoint,iCoord] = random.random() * 650;
        if (iCoord == 1):
            pointList[iPoint,iCoord] = random.random() * 450;
        if (iCoord == 2):
            pointList[iPoint,iCoord] = random.random() * 300;
writeVoxelizedImage ("someFleName.tif", pointList, 'parameter');
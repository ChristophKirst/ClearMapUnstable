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
    
    pointsToVTK(filename, x, y, z, data = data);
    
    
    
    # -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 15:17:01 2015

@author: kannanuv
"""
import numpy
import random


# vtkFileName is the file name to be written
# pointList is the list of x, y, z points to be written
# pointsLabel is the label for each point, it tell which region on brain the cell belongs to

def writeSparsePointsWithLabelVTKFile (vtkFileName, pointList, pointLabels):
    #swap x and y axis
    x = pointList[:,1];
    y = pointList[:,0];
    z = pointList[:,2];
    
    nPoint = x.size;
    
    #write VTK file
    vtkFilePointer = open(vtkFileName, 'w')
    vtkFilePointer.write('# vtk DataFile Version 2.0\n');
    vtkFilePointer.write('Unstructured Grid Example\n');
    vtkFilePointer.write('ASCII\n');
    vtkFilePointer.write('DATASET UNSTRUCTURED_GRID\n');
    vtkFilePointer.write("POINTS " + str(nPoint) + " float\n")
    for iPoint in range(0,nPoint):
        vtkFilePointer.write(str(x[iPoint]).format('%05.20f') + " " +  str(y[iPoint]).format('%05.20f') + " " + str(z[iPoint]).format('%05.20f') + "\n");    
    
    vtkFilePointer.write("CELLS " + str(nPoint) + " " + str(nPoint * 2) + "\n");
    for iPoint in range(0, nPoint):
        vtkFilePointer.write("1 " + str(iPoint) + "\n");
    vtkFilePointer.write("CELL_TYPES " + str(nPoint) + "\n");
    for iPoint in range(0, nPoint):
        vtkFilePointer.write("1 \n");
    #vtkFilePointer.write("\n");
    vtkFilePointer.write("POINT_DATA " + str(nPoint) + "\n");
    vtkFilePointer.write('SCALARS scalars float 1\n');
    vtkFilePointer.write("LOOKUP_TABLE default\n");
    for iLabel in pointLabels:
        vtkFilePointer.write(str(int(iLabel)) + " ");
        #vtkFilePointer.write("1 ")
    vtkFilePointer.write("\n");
    vtkFilePointer.close();
    
"""
    
    

    
    vtkFilePointer.close();
"""
# If the labels are not specified for the cells then all are assigned as one
def writeSparsePointsVTKFile (vtkFileName, pointList):
    x = pointList[:,1];
    nPoint = x.size;
    pointLabels = numpy.ones([nPoint,1], float)
    writeSparsePointsWithLabelVTKFile (vtkFileName, pointList, pointLabels)


# sample code of how to use the functions
vtkFileName = '/Users/kannanuv/Documents/workspace/code/assemblaSVN/idisco/iDISCO/Visualization/test.vtk'
nPoint = 100;
pointList = numpy.zeros([nPoint,3], float)
for iPoint in range(0, nPoint):
    for iCoord in range(0, 3):
        pointList[iPoint,iCoord] = random.random() * 10;
pointList
pointLabels = numpy.ones([nPoint,1], dtype=numpy.int)
writeSparsePointsVTKFile (vtkFileName, pointList)
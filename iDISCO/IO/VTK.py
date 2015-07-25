# -*- coding: utf-8 -*-
"""
Simple Interface to writing VTK files

Created on Thu Jun  4 14:37:06 2015

@author: ckirst / kannan
"""

#from evtk.hl import pointsToVTK 

import numpy;

#import iDISCO.IO.IO as io;

def writePoints(filename, points, labelimage = None):
    return "done"; 
#    #swap x and y axis
#    x = points[:,1];
#    y = points[:,0];
#    z = points[:,2];
##    
#    nPoint = x.size;
#    
#    if not labelimage is None:
#        if 
#        labelImage = io.readData(labelimage);
#        dsize =labelImage
#        for i in range(nPoint):
#            if x[i] > 0 and x[i] < dsize[0] and points[i,1] > 0 and points[i,1] < dsize[1] and points[i,2] > 0 and points[i,2] < dsize[2]:
#               imgc[points[i,0], points[i,1], points[i,2]] = 1;
#    self.plotOverlayLabel(image[:,:,0:dsize[2]], imgc, alpha = False);
#        pointLabels = labelImage, pointList(1,:), pointList(2,:) , pointList(3,:)
#    else
#        pointLabels = numpy.ones(nPoint);
#    end    
#    
##    
##    #write VTK file
#    vtkFilePointer = file(filename, 'w')
#    vtkFilePointer.write('# vtk DataFile Version 2.0\n');
#    vtkFilePointer.write('Unstructured Grid Example\n');
#    vtkFilePointer.write('ASCII\n');
#    vtkFilePointer.write('DATASET UNSTRUCTURED_GRID\n');
#    vtkFilePointer.write("POINTS " + str(nPoint) + " float\n")
#    for iPoint in range(0, nPoint):
#        vtkFilePointer.write(str(x[iPoint]).format('%05.20f') + " " +  str(y[iPoint]).format('%05.20f') + " " + str(z[iPoint]).format('%05.20f') + "\n");    
#    
#    vtkFilePointer.write("CELLS " + str(nPoint) + " " + str(nPoint * 2) + "\n");


#    for iPoint in range(0, nPoint):
#        vtkFilePointer.write("1 " + str(iPoint) + "\n");
#    vtkFilePointer.write("CELL_TYPES " + str(nPoint) + "\n");
#    for iPoint in range(0, nPoint):
#        vtkFilePointer.write("1 \n");
#    #vtkFilePointer.write("\n");
#    vtkFilePointer.write("POINT_DATA " + str(nPoint) + "\n");
#    vtkFilePointer.write('SCALARS scalars float 1\n');
#    vtkFilePointer.write("LOOKUP_TABLE default\n");
#    for iLabel in pointLabels:
#        vtkFilePointer.write(str(int(iLabel)) + " ");
#        #vtkFilePointer.write("1 ")
#    vtkFilePointer.write("\n");
#    vtkFilePointer.close();    
#
#
#
#
#if (numel(varargin) >= 2)
#    labelList = varargin{2};
#    labelIndices = ismember (pointLabels, labelList);
#    pointList = pointList (:, labelIndices);
#    pointLabels = pointLabels(labelIndices);
#end
#
#nPoint = size (pointList, 2);
#
#%%Swap x and y axis
#pointListx = pointList(2,:);
#pointListy = pointList(1,:);
#pointList(1,:) = pointListx;
#pointList(2,:) = pointListy;
#
#%% Write VTK File
#fp = fopen (vtkFileName, 'w+');
#fprintf (fp, '# vtk DataFile Version 2.0\n');
#fprintf (fp, 'Unstructured Grid Example\n');
#fprintf (fp, 'ASCII\n');
#fprintf (fp, 'DATASET UNSTRUCTURED_GRID\n');
#fprintf (fp, 'POINTS %d float\n', nPoint);
#fprintf (fp, '%05.20f %5.20f %5.20f\n', pointList);
#fprintf (fp, 'CELLS %d %d\n',nPoint, nPoint * 2);
#fprintf (fp, '%d %d\n', [ones(nPoint, 1) [0:nPoint-1]']');
#fprintf (fp, 'CELL_TYPES %d\n', nPoint);
#fprintf (fp, '%d\n', ones (nPoint, 1));
#fprintf (fp, 'POINT_DATA %d\n',nPoint);
#fprintf (fp, 'SCALARS scalars float 1\n');
#fprintf (fp, 'LOOKUP_TABLE default\n');
#fprintf (fp, '%f ', pointLabels);
#fclose (fp);
#
#
#
#
##def writePoints(filename, points, data = None):
##    x = points[:,0];
##    y = points[:,1];
##    z = points[:,2];
##    
##    pointsToVTK(filename, x, y, z, data = data);
#    
#  
#def readPoints(filename):
#    raise RuntimeError("Reading points from VTK file not implemented yet""")
#
#  
#  
#
# vtkFileName is the file name to be written
# pointList is the list of x, y, z points to be written
# pointsLabel is the label for each point, it tell which region on brain the cell belongs to
#
#def writeSparsePointsWithLabelVTKFile (vtkFileName, pointList, pointLabels):
#    #swap x and y axis
#    x = pointList[:,1];
#    y = pointList[:,0];
#    z = pointList[:,2];
#    
#    nPoint = x.size;
#    
#    #write VTK file
#    vtkFilePointer = open(vtkFileName, 'w')
#    vtkFilePointer.write('# vtk DataFile Version 2.0\n');
#    vtkFilePointer.write('Unstructured Grid Example\n');
#    vtkFilePointer.write('ASCII\n');
#    vtkFilePointer.write('DATASET UNSTRUCTURED_GRID\n');
#    vtkFilePointer.write("POINTS " + str(nPoint) + " float\n")
#    for iPoint in range(0,nPoint):
#        vtkFilePointer.write(str(x[iPoint]).format('%05.20f') + " " +  str(y[iPoint]).format('%05.20f') + " " + str(z[iPoint]).format('%05.20f') + "\n");    
#    
#    vtkFilePointer.write("CELLS " + str(nPoint) + " " + str(nPoint * 2) + "\n");
#    for iPoint in range(0, nPoint):
#        vtkFilePointer.write("1 " + str(iPoint) + "\n");
#    vtkFilePointer.write("CELL_TYPES " + str(nPoint) + "\n");
#    for iPoint in range(0, nPoint):
#        vtkFilePointer.write("1 \n");
#    #vtkFilePointer.write("\n");
#    vtkFilePointer.write("POINT_DATA " + str(nPoint) + "\n");
#    vtkFilePointer.write('SCALARS scalars float 1\n');
#    vtkFilePointer.write("LOOKUP_TABLE default\n");
#    for iLabel in pointLabels:
#        vtkFilePointer.write(str(int(iLabel)) + " ");
#        #vtkFilePointer.write("1 ")
#    vtkFilePointer.write("\n");
#    vtkFilePointer.close();
#    
#"""
#    
#    
#
#    
#    vtkFilePointer.close();
#"""
## If the labels are not specified for the cells then all are assigned as one
#def writeSparsePointsVTKFile (vtkFileName, pointList):
#    x = pointList[:,1];
#    nPoint = x.size;
#    pointLabels = numpy.ones([nPoint,1], float)
#    writeSparsePointsWithLabelVTKFile (vtkFileName, pointList, pointLabels)
#
#
## sample code of how to use the functions
#vtkFileName = '/Users/kannanuv/Documents/workspace/code/assemblaSVN/idisco/iDISCO/Visualization/test.vtk'
#nPoint = 100;
#pointList = numpy.zeros([nPoint,3], float)
#for iPoint in range(0, nPoint):
#    for iCoord in range(0, 3):
#        pointList[iPoint,iCoord] = random.random() * 10;
#pointList
#pointLabels = numpy.ones([nPoint,1], dtype=numpy.int)
#writeSparsePointsVTKFile (vtkFileName, pointList)
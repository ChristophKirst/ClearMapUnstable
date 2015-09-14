# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 17:43:39 2015

@author: mtllab
"""

import iDISCO.IO.IO as io
import numpy

from scipy.ndimage.measurements import label

import os, numpy, math
import matplotlib.pyplot as pp
import iDISCO.Visualization.Plot as plt;
import iDISCO.Settings as settings
from iDISCO.Alignment.Resampling import resampleData;
from iDISCO.Alignment.Elastix import alignData, transformPoints
from iDISCO.ImageProcessing.CellDetection import detectCells
from iDISCO.Alignment.Resampling import resamplePoints, resamplePointsInverse
from iDISCO.Analysis.Label import countPointsInRegions
from iDISCO.Analysis.Statistics import thresholdPoints
from iDISCO.Utils.ParameterTools import joinParameter
import iDISCO.Analysis.Label as lbl
from iDISCO.Analysis.Voxelization import voxelize

from scipy.ndimage.morphology import binary_dilation


annotationFiles = ['/home/mtllab/Documents/countingaccuracy/cortex/nico.nrrd',
              '/home/mtllab/Documents/countingaccuracy/cortex/eliza.nrrd'];

#annotationData = [io.readData(l, y = (50, 151)) for l in annotationFiles];
annotationData = [io.readData(l) for l in annotationFiles];


annotationLabel = [label(l) for l in annotationData];
annotationNLabel =  [l[1] for l in annotationLabel];
annotationLabel = [l[0] for l in annotationLabel];


correspondence = [];
multipleOverlap1 = [];
noOverlap1 = [];

for l in range(1, annotationNLabel[0] + 1):
  #  print l
    dataL1 = annotationLabel[0] == l;
    
    ll = annotationLabel[1][dataL1];
    ll = numpy.hstack((ll, numpy.array(0)));
    
    if numpy.any(ll):
        q, c = numpy.unique(ll, return_counts = True);
        for i in range(q.shape[0]):
            p = q[i];
            if p > 0:
                correspondence.append((l,p, c[i]));
                
        if q.shape[0] > 2:
            multipleOverlap1.append((l, q[1:]));   
    else:
        noOverlap1.append(l);
        
        
c = numpy.array( [cc[1] for cc in correspondence] );  
noOverlap2 = numpy.setdiff1d(range(1, annotationNLabel[1] + 1), c);


c = numpy.array( [cc[1] for cc in correspondence] );
a, b = numpy.unique(c, return_counts = True);
mo2 = numpy.sum(b- 1)


print  annotationNLabel  
print "#pairs:%d,  multiple overlap1:%d,  no overlap1 %d, multiple overlap2: %d, no overlap2 %d" % (len(correspondence), len(multipleOverlap1), len(noOverlap1), mo2, noOverlap2.shape[0]);






from scipy.ndimage.morphology import binary_dilation

#points = io.readPoints('/home/mtllab/Documents/countingaccuracy/cells-40.npy');

#annotationCode = io.readData('/home/mtllab/Documents/countingaccuracy/cortex/cellsnico_dog446-40.tif');
annotationCode = io.readData('/home/mtllab/Documents/countingaccuracy/cortex/cells_unfiltered700.tif');


annotationCodeD = binary_dilation(annotationCode, numpy.ones((3,3,3), dtype = bool));
annotationCodeD = annotationCodeD.astype('int32');

io.writeData('/home/mtllab/Desktop/testD.nrrd', annotationCodeD)


annotationLabelCode = label(annotationCodeD);
annotationLabelCodeN = annotationLabelCode[1];
annotationLabelCode = annotationLabelCode[0];


correspondence = [];
multipleOverlap1 = [];
noOverlap1 = [];

for l in range(1, annotationNLabel[0] + 1):
  #  print l
    dataL1 = annotationLabel[0] == l;
    
    ll = annotationLabelCode[dataL1];
    ll = numpy.hstack((ll, numpy.array(0)));
    
    if numpy.any(ll):
        q, c = numpy.unique(ll, return_counts = True);
        for i in range(q.shape[0]):
            p = q[i];
            if p > 0:
                correspondence.append((l,p, c[i]));
                
        if q.shape[0] > 2:
            multipleOverlap1.append((l, q[1:]));   
    else:
        noOverlap1.append(l);
        
        
c = numpy.array( [cc[1] for cc in correspondence] );  
noOverlap2 = numpy.setdiff1d(range(1, annotationLabelCodeN + 1), c);


c = numpy.array( [cc[1] for cc in correspondence] );
a, b = numpy.unique(c, return_counts = True);
mo2 = numpy.sum(b- 1)


print  annotationNLabel[0], annotationLabelCodeN
print "#pairs:%d,  multiple overlap1:%d,  no overlap1 %d, multiple overlap2: %d, no overlap2 %d" % (len(correspondence), len(multipleOverlap1), len(noOverlap1), mo2, noOverlap2.shape[0]);






correspondence = [];
multipleOverlap1 = [];
noOverlap1 = [];

for l in range(1, annotationNLabel[1] + 1):
   # print l
    dataL1 = annotationLabel[1] == l;
    
    ll = annotationLabelCode[dataL1];
    ll = numpy.hstack((ll, numpy.array(0)));
    
    if numpy.any(ll):
        q, c = numpy.unique(ll, return_counts = True);
        for i in range(q.shape[0]):
            p = q[i];
            if p > 0:
                correspondence.append((l,p, c[i]));
                
        if q.shape[0] > 2:
            multipleOverlap1.append((l, q[1:]));   
    else:
        noOverlap1.append(l);
        
        
c = numpy.array( [cc[1] for cc in correspondence] );  
noOverlap2 = numpy.setdiff1d(range(1, annotationLabelCodeN + 1), c);


c = numpy.array( [cc[1] for cc in correspondence] );
a, b = numpy.unique(c, return_counts = True);
mo2 = numpy.sum(b- 1)


print  annotationNLabel[0], annotationLabelCodeN
print "#pairs:%d,  multiple overlap1:%d,  no overlap1 %d, multiple overlap2: %d, no overlap2 %d" % (len(correspondence), len(multipleOverlap1), len(noOverlap1), mo2, noOverlap2.shape[0]);



#######################################################



BaseDirectory = '/home/mtllab/Documents/countingaccuracy/';
   
points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells_dog-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities_dog-allpoints.npy')));

annotationFiles = ['/home/mtllab/Documents/countingaccuracy/nico.nrrd']
            #  '/home/mtllab/Documents/countingaccuracy/nico.nrrd'];

#annotationData = [io.readData(l, y = (50, 151)) for l in annotationFiles];
annotationData = [io.readData(l) for l in annotationFiles];


annotationLabel = [label(l) for l in annotationData];
annotationNLabel =  [l[1] for l in annotationLabel];
annotationLabel = [l[0] for l in annotationLabel];


correspondence = [];
multipleOverlap1 = [];
noOverlap1 = [];

dataSource =  '/home/mtllab/Documents/countingaccuracy/countingaccuracy.tif'

## Parameter to calculate density voxelization
VoxelizationParameter = {
    #Method to voxelize
    "voxelizationMethod" : 'Pixel', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized in pixels
    "voxelizationSize" : (1,1,1),  

    # Voxelization weigths (e/g intensities)
    "voxelizationWeights" : None
    };
check = [0,0,0,0,0];

for t in range(0,150):
    correspondence = [];
    multipleOverlap1 = [];
    noOverlap1 = [];    
    pointT, intensitiesT = thresholdPoints(points, intensities, threshold = (t, 65000), row = (1,0));
    annotationCode = voxelize(pointT, dataSource, **VoxelizationParameter);
    annotationCodeD = binary_dilation(annotationCode, numpy.ones((3,3,3), dtype = bool));
    annotationCodeD = annotationCodeD.astype('int32');
    annotationLabelCode = label(annotationCodeD);    
    annotationLabelCodeN = annotationLabelCode[1];
    annotationLabelCode = annotationLabelCode[0];
    for l in range(1, annotationNLabel[0] + 1):
        dataL1 = annotationLabel[0] == l;
        ll = annotationLabelCode[dataL1];
        ll = numpy.hstack((ll, numpy.array(0)));
    
        if numpy.any(ll):
            q, c = numpy.unique(ll, return_counts = True);
            for i in range(q.shape[0]):
                p = q[i];
                if p > 0:
                    correspondence.append((l,p, c[i]));
                
                if q.shape[0] > 2:
                    multipleOverlap1.append((l, q[1:]));   
                else:
                    noOverlap1.append(l);
        
        
        c = numpy.array( [cc[1] for cc in correspondence] );  
        noOverlap2 = numpy.setdiff1d(range(1, annotationLabelCodeN + 1), c);


        c = numpy.array( [cc[1] for cc in correspondence] );
        a, b = numpy.unique(c, return_counts = True);
        mo2 = numpy.sum(b- 1)
    print t
    print  annotationNLabel[0], annotationLabelCodeN
    print "#pairs:%d,  multiple overlap1:%d,  no overlap1 %d, multiple overlap2: %d, no overlap2 %d" % (len(correspondence), len(multipleOverlap1), len(noOverlap1), mo2, noOverlap2.shape[0]);
    check =  numpy.vstack([check,[t,annotationNLabel[0],annotationLabelCodeN,len(correspondence),noOverlap2.shape[0]]]);


io.writePoints('/home/mtllab/Documents/countingaccuracy/nicocheck_dog446.csv', check);



#######################################################



BaseDirectory = '/home/mtllab/Documents/countingaccuracy/';
   
points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells_dog-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities_dog-allpoints.npy')));

annotationFiles = ['/home/mtllab/Documents/countingaccuracy/nico.nrrd']
            #  '/home/mtllab/Documents/countingaccuracy/nico.nrrd'];

#annotationData = [io.readData(l, y = (50, 151)) for l in annotationFiles];
annotationData = [io.readData(l) for l in annotationFiles];


annotationLabel = [label(l) for l in annotationData];
annotationNLabel =  [l[1] for l in annotationLabel];
annotationLabel = [l[0] for l in annotationLabel];


correspondence = [];
multipleOverlap1 = [];
noOverlap1 = [];

dataSource =  '/home/mtllab/Documents/countingaccuracy/countingaccuracy.tif'

## Parameter to calculate density voxelization
VoxelizationParameter = {
    #Method to voxelize
    "voxelizationMethod" : 'Pixel', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized in pixels
    "voxelizationSize" : (1,1,1),  

    # Voxelization weigths (e/g intensities)
    "voxelizationWeights" : None
    };
check = [0,0,0,0,0];

for t in range(0,60):
    correspondence = [];
    multipleOverlap1 = [];
    noOverlap1 = [];    
    pointT, intensitiesT = thresholdPoints(points, intensities, threshold = (t, 65000), row = (1,0));
    annotationCode = voxelize(pointT, dataSource, **VoxelizationParameter);
    annotationCodeD = binary_dilation(annotationCode, numpy.ones((3,3,3), dtype = bool));
    annotationCodeD = annotationCodeD.astype('int32');
    annotationLabelCode = label(annotationCodeD);    
    annotationLabelCodeN = annotationLabelCode[1];
    annotationLabelCode = annotationLabelCode[0];
    for l in range(1, annotationNLabel[0] + 1):
        dataL1 = annotationLabel[0] == l;
        ll = annotationLabelCode[dataL1];
        ll = numpy.hstack((ll, numpy.array(0)));
    
        if numpy.any(ll):
            q, c = numpy.unique(ll, return_counts = True);
            for i in range(q.shape[0]):
                p = q[i];
                if p > 0:
                    correspondence.append((l,p, c[i]));
                
                if q.shape[0] > 2:
                    multipleOverlap1.append((l, q[1:]));   
                else:
                    noOverlap1.append(l);
        
        
        c = numpy.array( [cc[1] for cc in correspondence] );  
        noOverlap2 = numpy.setdiff1d(range(1, annotationLabelCodeN + 1), c);


        c = numpy.array( [cc[1] for cc in correspondence] );
        a, b = numpy.unique(c, return_counts = True);
        mo2 = numpy.sum(b- 1)
    print t
    print  annotationNLabel[0], annotationLabelCodeN
    print "#pairs:%d,  multiple overlap1:%d,  no overlap1 %d, multiple overlap2: %d, no overlap2 %d" % (len(correspondence), len(multipleOverlap1), len(noOverlap1), mo2, noOverlap2.shape[0]);
    check =  numpy.vstack([check,[t,annotationNLabel[0],annotationLabelCodeN,len(correspondence),noOverlap2.shape[0]]]);


io.writePoints('/home/mtllab/Documents/countingaccuracy/elizacheck_dog.csv', check);


####################################################################################
####################################################################################
#######################################################




######################### Data parameters

#BaseDirectory = '/home/mtllab/Documents/countingaccuracy/cortex';

BaseDirectory = '/home/mtllab/Documents/countingaccuracy/cortex';
  
  
#cFosFile = '/home/mtllab/Documents/countingaccuracy/cortex/data.tif';
cFosFile =  '/home/mtllab/Documents/countingaccuracy/cortex/data.tif'

cFosFileRange = {'x' : all, 'y' : all, 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (815,1000), 'y' : (1078,1271), 'z' : (667,742)};


#Resolution of the Data (in um / pixel)
OriginalResolution = (4.0625, 4.0625, 3);


######################### Cell Detection Parameters using DoG

SpotDetectionParameter = {
    # background correctoin: None or (x,y) which is size of disk for gray scale opening
    "backgroundSize" : (7,7),
    
    # spot Detection via Difference of Gaussians (DoG) filter: (x,y,z) size
    "dogSize" : None,
    
    # h of h-max transform
    "hMax" : None,
    
    # intensity detection   
    "intensityMethod"  : 'Max',  #None -> intensity of pixel of center, alternatively string of numpy array method that returns a single number
    "intensitySize"    : (3,3,3),  # size of box in (x,y,z) to include in intensity determination
    
    # threshold for min intensity at center to be counted as cell, for saving ('None' will save everything )
    "threshold" : None,
      
    # write cell mask to disk (to check cell detection accuracy), if not None
    #"cellMaskFile" : os.path.join(BaseDirectory, 'cell_mask/cell_mask_Z\d{4}.ome.tif')
    "cellMaskFile" : None
    
    #some debug / quality check output
    #"verbose" : True,
    #"processMethod" : "sequential"  #  plotting during image processing only in sequential mode !
    };



############################ Config parameters
VoxelizationParameter = {
    #Method to voxelize
    "voxelizationMethod" : 'Pixel', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized in pixels
    "voxelizationSize" : (1,1,1),  

    # Voxelization weigths (e/g intensities)
    "voxelizationWeights" : None
    };
#Processes to use for Resampling
ResamplingParameter = { "processes": 12 };


#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 10,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "chunkOverlap" : 30,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : False,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
    "chunkOptimizationSize" : all
    };


# result files for cell coordinates (csv, vtk or ims)
SpotDetectionParameter["source"] = cFosFile;
SpotDetectionParameter = joinParameter(SpotDetectionParameter, cFosFileRange)

SpotDetectionParameter["sink"] = (None, None);

ImageProcessingParameter = joinParameter(StackProcessingParameter, SpotDetectionParameter)


annotationFiles = ['/home/mtllab/Documents/countingaccuracy/cortex/nico.nrrd']
            #  '/home/mtllab/Documents/countingaccuracy/nico.nrrd'];

#annotationData = [io.readData(l, y = (50, 151)) for l in annotationFiles];
annotationData = [io.readData(l) for l in annotationFiles];


annotationLabel = [label(l) for l in annotationData];
annotationNLabel =  [l[1] for l in annotationLabel];
annotationLabel = [l[0] for l in annotationLabel];


correspondence = [];
multipleOverlap1 = [];
noOverlap1 = [];

dataSource =  '/home/mtllab/Documents/countingaccuracy/cortex/data.tif'

check = [0,0,0,0,0,0,0,0,0];


for dogxy in range(1,10):
  #  for dogz in range(1,10):
        SpotDetectionParameter["dogSize"] = (dogxy,dogxy,2*dogxy);
        ImageProcessingParameter = joinParameter(StackProcessingParameter, SpotDetectionParameter)
        points, intensities = detectCells(**ImageProcessingParameter);
        
        for t in range(0,70):
            correspondence = [];
            multipleOverlap1 = [];
            noOverlap1 = [];    
            pointT, intensitiesT = thresholdPoints(points, intensities, threshold = (t, 65000), row = (1,0));
            annotationCode = voxelize(pointT, dataSource, **VoxelizationParameter);
            annotationCodeD = binary_dilation(annotationCode, numpy.ones((3,3,3), dtype = bool));
            annotationCodeD = annotationCodeD.astype('int32');
            annotationLabelCode = label(annotationCodeD);    
            annotationLabelCodeN = annotationLabelCode[1];
            annotationLabelCode = annotationLabelCode[0];
            for l in range(1, annotationNLabel[0] + 1):
                dataL1 = annotationLabel[0] == l;
                ll = annotationLabelCode[dataL1];
                ll = numpy.hstack((ll, numpy.array(0)));
            
                if numpy.any(ll):
                    q, c = numpy.unique(ll, return_counts = True);
                    for i in range(q.shape[0]):
                        p = q[i];
                        if p > 0:
                            correspondence.append((l,p, c[i]));
                        
                        if q.shape[0] > 2:
                            multipleOverlap1.append((l, q[1:]));   
                        else:
                            noOverlap1.append(l);
                
                
                c = numpy.array( [cc[1] for cc in correspondence] );  
                noOverlap2 = numpy.setdiff1d(range(1, annotationLabelCodeN + 1), c);
        
        
                c = numpy.array( [cc[1] for cc in correspondence] );
                a, b = numpy.unique(c, return_counts = True);
                mo2 = numpy.sum(b- 1)
            print t
            print  annotationNLabel[0], annotationLabelCodeN
            print "#pairs:%d,  multiple overlap1:%d,  no overlap1 %d, multiple overlap2: %d, no overlap2 %d" % (len(correspondence), len(multipleOverlap1), len(noOverlap1), mo2, noOverlap2.shape[0]);
            check =  numpy.vstack([check,[dogxy, dogz,t,annotationNLabel[0],annotationLabelCodeN,len(correspondence),noOverlap2.shape[0],(len(correspondence)/float(annotationNLabel[0])),(noOverlap2.shape[0]/float(annotationNLabel[0]))]]);


io.writePoints('/home/mtllab/Documents/countingaccuracy/cortex/nicocheck_dogz.csv', check);


##########################################################################


######################### Data parameters

#BaseDirectory = '/home/mtllab/Documents/countingaccuracy/cortex';

BaseDirectory = '/home/mtllab/Documents/countingaccuracy/cortex';
  
  
#cFosFile = '/home/mtllab/Documents/countingaccuracy/cortex/data.tif';
cFosFile =  '/home/mtllab/Documents/countingaccuracy/cortex/data.tif'

cFosFileRange = {'x' : all, 'y' : all, 'z' : all};#cFosFileRange = {'x' : all, 'y' : (180, 2560), 'z' : all};
#cFosFileRange = {'x' : (815,1000), 'y' : (1078,1271), 'z' : (667,742)};


#Resolution of the Data (in um / pixel)
OriginalResolution = (4.0625, 4.0625, 3);


######################### Cell Detection Parameters using DoG

SpotDetectionParameter = {
    # background correctoin: None or (x,y) which is size of disk for gray scale opening
    "backgroundSize" : (15,15),
    
    # spot Detection via Difference of Gaussians (DoG) filter: (x,y,z) size
    "dogSize" : None,
    
    # h of h-max transform
    "hMax" : None,
    
    # intensity detection   
    "intensityMethod"  : None,  #None -> intensity of pixel of center, alternatively string of numpy array method that returns a single number
    "intensitySize"    : (3,3,3),  # size of box in (x,y,z) to include in intensity determination
    
    # threshold for min intensity at center to be counted as cell, for saving ('None' will save everything )
    "threshold" : None,
      
    # write cell mask to disk (to check cell detection accuracy), if not None
    #"cellMaskFile" : os.path.join(BaseDirectory, 'cell_mask/cell_mask_Z\d{4}.ome.tif')
    "cellMaskFile" : None
    
    #some debug / quality check output
    #"verbose" : True,
    #"processMethod" : "sequential"  #  plotting during image processing only in sequential mode !
    };



############################ Config parameters
VoxelizationParameter = {
    #Method to voxelize
    "voxelizationMethod" : 'Pixel', # Spherical,'Rectangular, Gaussian'
       
    # Define bounds of the volume to be voxelized in pixels
    "voxelizationSize" : (1,1,1),  

    # Voxelization weigths (e/g intensities)
    "voxelizationWeights" : None
    };
#Processes to use for Resampling
ResamplingParameter = { "processes": 12 };


#Stack Processing Parameter for cell detection
StackProcessingParameter = {
    #max number of parallel processes
    "processes" : 10,
   
    #chunk sizes
    "chunkSizeMax" : 100,
    "chunkSizeMin" : 30,
    "chunkOverlap" : 30,

    #optimize chunk size and number to number of processes
    "chunkOptimization" : False,
    
    #increase chunk size for optimizaition (True, False or all = automatic)
    "chunkOptimizationSize" : all
    };


# result files for cell coordinates (csv, vtk or ims)
SpotDetectionParameter["source"] = cFosFile;
SpotDetectionParameter = joinParameter(SpotDetectionParameter, cFosFileRange)

SpotDetectionParameter["sink"] = (None, None);

ImageProcessingParameter = joinParameter(StackProcessingParameter, SpotDetectionParameter)


annotationFiles = ['/home/mtllab/Documents/countingaccuracy/cortex/nico.nrrd']
            #  '/home/mtllab/Documents/countingaccuracy/nico.nrrd'];

#annotationData = [io.readData(l, y = (50, 151)) for l in annotationFiles];
annotationData = [io.readData(l) for l in annotationFiles];


annotationLabel = [label(l) for l in annotationData];
annotationNLabel =  [l[1] for l in annotationLabel];
annotationLabel = [l[0] for l in annotationLabel];


correspondence = [];
multipleOverlap1 = [];
noOverlap1 = [];

dataSource =  '/home/mtllab/Documents/countingaccuracy/cortex/data.tif'

check = [0,0,0,0,0,0,0,0,0,0];

dogxy = 0;
dogz = 0;

for bgd in range(3,40,2):
  #  for dogz in range(3,10):
     #   for dogxy in range(3,10):        
            SpotDetectionParameter['backgroundSize'] = (bgd,bgd);
         #   SpotDetectionParameter['dogSize'] = (dogxy,dogxy,dogz);
            ImageProcessingParameter = joinParameter(StackProcessingParameter, SpotDetectionParameter)
            points, intensities = detectCells(**ImageProcessingParameter);
            
            for t in range(0,10000,100):
                correspondence = [];
                multipleOverlap1 = [];
                noOverlap1 = [];    
                pointT, intensitiesT = thresholdPoints(points, intensities, threshold = (t, 65000), row = (1,0));
                annotationCode = voxelize(pointT, dataSource, **VoxelizationParameter);
                annotationCodeD = binary_dilation(annotationCode, numpy.ones((3,3,3), dtype = bool));
                annotationCodeD = annotationCodeD.astype('int32');
                annotationLabelCode = label(annotationCodeD);    
                annotationLabelCodeN = annotationLabelCode[1];
                annotationLabelCode = annotationLabelCode[0];
                for l in range(1, annotationNLabel[0] + 1):
                    dataL1 = annotationLabel[0] == l;
                    ll = annotationLabelCode[dataL1];
                    ll = numpy.hstack((ll, numpy.array(0)));
                
                    if numpy.any(ll):
                        q, c = numpy.unique(ll, return_counts = True);
                        for i in range(q.shape[0]):
                            p = q[i];
                            if p > 0:
                                correspondence.append((l,p, c[i]));
                            
                            if q.shape[0] > 2:
                                multipleOverlap1.append((l, q[1:]));   
                            else:
                                noOverlap1.append(l);
                    
                    
                    c = numpy.array( [cc[1] for cc in correspondence] );  
                    noOverlap2 = numpy.setdiff1d(range(1, annotationLabelCodeN + 1), c);
            
            
                    c = numpy.array( [cc[1] for cc in correspondence] );
                    a, b = numpy.unique(c, return_counts = True);
                    mo2 = numpy.sum(b- 1)
                print t
                print  annotationNLabel[0], annotationLabelCodeN
                print "#pairs:%d,  multiple overlap1:%d,  no overlap1 %d, multiple overlap2: %d, no overlap2 %d" % (len(correspondence), len(multipleOverlap1), len(noOverlap1), mo2, noOverlap2.shape[0]);
                check =  numpy.vstack([check,[bgd,dogxy,dogz,t,annotationNLabel[0],annotationLabelCodeN,len(correspondence),noOverlap2.shape[0],(len(correspondence)/float(annotationNLabel[0])),(noOverlap2.shape[0]/float(annotationNLabel[0]))]]);
    

io.writePoints('/home/mtllab/Documents/countingaccuracy/cortex/nicocheck_parameters_nofilters.csv', check);



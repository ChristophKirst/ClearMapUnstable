# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 17:43:39 2015

@author: mtllab
"""

import iDISCO.IO.IO as io
import numpy

import iDISCO.Analysis.Voxelization as vox

from scipy.ndimage.measurements import label



annotationFiles = ['/home/mtllab/Documents/countingaccuracy/eliza.tif',
              '/home/mtllab/Documents/countingaccuracy/nico.nrrd'];

#annotationData = [io.readData(l, y = (50, 151)) for l in annotationFiles];
annotationData = [io.readData(l) for l in annotationFiles];


annotationLabel = [label(l) for l in annotationData];
annotationNLabel =  [l[1] for l in annotationLabel];
annotationLabel = [l[0] for l in annotationLabel];


correspondence = [];
multipleOverlap1 = [];
noOverlap1 = [];

for l in range(1, annotationNLabel[0] + 1):
    print l
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

annotationCode = io.readData('/home/mtllab/Documents/countingaccuracy/countingaccuracy-30.tif');

annotationCode = io.readData('/home/mtllab/Documents/countingaccuracy/countingaccuracy-2200_fast.tif');


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
    print l
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
    print l
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









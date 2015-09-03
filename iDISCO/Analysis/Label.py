# -*- coding: utf-8 -*-
"""
Label and Annotation info from Allen Brain Atlas (v2)

Created on Thu Aug  6 19:25:15 2015

@author: ckirst
"""

import numpy
import iDISCO.IO.IO as io

import sys
self = sys.modules[__name__];

import os
import csv

from iDISCO.Settings import IDISCOPath

import collections

LabelRecord = collections.namedtuple('LabelRecord', 'id, name, acronym, color, parent');

DefaultLabeledImageFile = os.path.join(IDISCOPath, 'Test/Data/Annotation/annotation_25_right.tif');
DefaultAnnotationFile   = os.path.join(IDISCOPath, 'Data/ARA2_annotation_info.csv');



def labelToInt(txt):
    if txt == '':
        return -1;
    else:
        return int(txt);

class LabelInfo(object):
    ids = None;
    names = None;
    acronyms = None;
    colors = None;
    parents = None;
    levels = None;
    
    def __init__(slf, annotationFile = DefaultAnnotationFile):
        slf.initialize(annotationFile = annotationFile);
    
    def initialize(slf, annotationFile = DefaultAnnotationFile):
        
        with open(annotationFile) as dfile:
            reader = csv.reader(dfile);
            #skip header
            reader.next();
            labels = [LabelRecord._make((int(row[0]), row[1], row[2], [int(row[3]), int(row[4]), int(row[5])], labelToInt(row[7]))) for row in reader];
            
            dfile.close();
        
        slf.ids = [x.id for x in labels];
        slf.names = { x.id : x.name for x in labels};
        slf.acronyms = { x.id : x.acronym for x in labels};
        slf.colors = { x.id : x.color for x in labels};
        slf.parents = { x.id : x.parent for x in labels};
        
        #calculate level:
        slf.levels = [0 for x in labels];
        for i in range(len(slf.ids)):
            p = slf.ids[i];
            slf.levels[i] = 0;
            p = slf.parents[p];
            while p >= 0:
                slf.levels[i] += 1;
                p = slf.parents[p];
        slf.levels = {slf.ids[i] : slf.levels[i] for i in range(len(slf.levels))};
            
    
    def name(slf, iid):
        return slf.names[iid];
        
    def acronym(slf, iid):
        return slf.acronyms[iid];
        
    def color(slf, iid):
        return slf.colors[iid];
        
    def parent(slf, iid):
        return slf.parents[iid];
        
    def level(slf, iid):
        return slf.levels[iid];
    
    def toLabelAtLevel(slf, iid, level):
        i = iid;
        if not i in slf.ids:
            return i;
        
        #print i, levels[i]
        while slf.levels[i] > level:
            i = slf.parents[i];
            #print i, levels[i]
        return i;
            

Label = LabelInfo();

def initialize(annotationFile = DefaultAnnotationFile):
    Label.initialize(annotationFile);


def labelAtLevel(label, level):
    global Label;
    
    if level is None:
        return label;
    else:
        if isinstance(label, numpy.ndarray):
            return [Label.toLabelAtLevel(x, level) for x in label];
        else:
            return Label.toLabelAtLevel(label, level);


def labelPoints(points, labeledImage = DefaultLabeledImageFile, level = None):
    
    #points are (y,x,z) -> which is also the way the labeled image is read in
    #x = points[:,1];
    #y = points[:,0];
    #z = points[:,2];

    x = points[:,0];
    y = points[:,1];
    z = points[:,2]; 
    
    nPoint = x.size;    
    
    pointLabels = numpy.zeros(nPoint, 'int32');
    
    labelImage = io.readData(labeledImage);
                  
    dsize = labelImage.shape;
    for i in range(nPoint):
        #if y[i] >= 0 and y[i] < dsize[0] and x[i] >= 0 and x[i] < dsize[1] and z[i] >= 0 and z[i] < dsize[2]:
        #     pointLabels[i] = labelImage[y[i], x[i], z[i]];
        if x[i] >= 0 and x[i] < dsize[0] and y[i] >= 0 and y[i] < dsize[1] and z[i] >= 0 and z[i] < dsize[2]:
             pointLabels[i] = labelImage[x[i], y[i], z[i]];
    
    pointLabels = self.labelAtLevel(pointLabels, level);
    
    return pointLabels;


 
def countPointsInRegions(points, labeledImage = DefaultLabeledImageFile, intensities = None, intensityRow = 0, level= None, allIds = False, sort = True, returnIds = True, returnCounts = False):
    global Label;
    
    points = io.readPoints(points);
    intensities = io.readPoints(intensities);
    pointLabels = self.labelPoints(points, labeledImage, level = level); 
    
    if intensities is None:
        ll, cc = numpy.unique(pointLabels, return_counts = True);
        cci = None;
    else:
        if intensities.ndim > 1:
            intensities = intensities[:,intensityRow];
   
        ll, ii, cc = numpy.unique(pointLabels, return_counts = True, return_inverse = True);
        cci = numpy.zeros(ll.shape);
        for i in range(ii.shape[0]):
             cci[ii[i]] += intensities[i];
    
    if allIds:
        lla = numpy.setdiff1d(Label.ids, ll);
        ll  = numpy.hstack((ll, lla));
        cc  = numpy.hstack((cc, numpy.zeros(lla.shape, dtype = cc.dtype)));
        if not cci is None:
            cci = numpy.hstack((cci, numpy.zeros(lla.shape, dtype = cc.dtype)));
        
    
    #cc = numpy.vstack((ll,cc)).T;
    if sort:
        ii = numpy.argsort(ll);
        cc = cc[ii];
        ll = ll[ii];
        if not cci is None:
            cci = cci[ii];

    if returnIds:
        if cci is None:
            return ll, cc
        else:
            if returnCounts:
                return ll, cc, cci;
            else:
                return ll, cci
    else:
        if cci is None:
            return cc;
        else:
            if returnCounts:
                return cc, cci;
            else:
                return cci;

    
def labelToName(label):
    global Label;
    return [Label.name(x) for x in label];
    
def labelToAcronym(label):
    global Label;
    return [Label.acronym(x) for x in label];

def labelToColor(label):
    global Label;
    return [Label.color(x) for x in label];



#def writeCountsToLabeledImage(counts, labelimage = defaultLabeledImage(), level = None):
#    #read the labeled image
#    labelimage = io.readData(labelimage);
 

   

def writePAL(filename, cols):
    with open(filename, 'w') as f:
        for c in cols:
            f.write(('%3.3f' % c[0]).rjust(10) +  ('%3.3f' % c[1]).rjust(11) + ('%3.3f' % c[2]).rjust(11));
            f.write('\n');
        f.close();

def writeLUT(filename, cols):
    with open(filename, 'w') as f:
        for c in cols:
            f.write(('%d' % c[0]).rjust(4) +  ('%d' % c[1]).rjust(5) + ('%d' % c[2]).rjust(5));
            f.write('\n');
        f.close();

      
def makeColorPalette(filename = None):
    """Creates a pal file for imaris based on label colors"""
    
    global Label; 
    maxlabel = max(Label.ids);
    colarray = numpy.zeros((maxlabel, 3));
    for i in Label.ids:
        colarray[i-1,:] = Label.color(i);
    
    if filename is None:
        return colarray;
    else:
        fext = io.fileExtension(filename);
        if fext == 'pal':
            writePAL(filename, colarray);
        elif fext == 'lut':
            writeLUT(filename, colarray);
        else:
            raise RuntimeError('color pallete format: %s not supported' % fext);
        
        return colarray;
            

def makeColorAnnotations(filename, labeledImage = None):
    if labeledImage is None:
        labeledImage = DefaultLabeledImageFile;
    
    li = io.readData(labeledImage);
    dsize = li.shape;
    
    lr = numpy.zeros(dsize, dtype = numpy.uint8);
    lg = lr.copy();
    lb = lr.copy();
    
    global Label; 
    maxlabel = max(Label.ids);
    colarray = numpy.zeros((maxlabel, 3));
    for i in Label.ids:
        colarray[i-1,:] = Label.color(i);
    
    for i in Label.ids:
        ll = li == i;
        lr[ll] = colarray[i-1,0];
        lg[ll] = colarray[i-1,1];
        lb[ll] = colarray[i-1,2];
    
    io.writeData(filename + "_r.tif", lr);
    io.writeData(filename + "_g.tif", lg);  
    io.writeData(filename + "_b.tif", lb);
    
    return (lr,lg,lb);


def test():
    """Test Label module"""

    import iDISCO.Analysis.Label as self
    import numpy
    points = numpy.array([[162, 200, 138], [246, 486, 138], [246, 486, 138]]);
    intensities = numpy.array([1,2,3]);
    
    ll, cc = self.countPointsInRegions(points, intensities, withIds = True);
    


if __name__ == '__main__':
    test();
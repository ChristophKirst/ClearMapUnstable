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

from iDISCO.Parameter import iDISCOPath


def toInt(txt):
    if txt == '':
        return -1;
    else:
        return int(txt);
        
import collections
LabelRecord = collections.namedtuple('LabelRecord', 'id, name, acronym, color, parent');

class LabelInfo(object):
    ids = None;
    names = None;
    acronyms = None;
    colors = None;
    parents = None;
    levels = None;
    
    def initialize(slf):
        datafile = os.path.join(iDISCOPath(), 'Data/ARA2_annotation_info.csv');
        
        with open(datafile) as dfile:
            reader = csv.reader(dfile);
            #skip header
            reader.next();
            labels = [LabelRecord._make((int(row[0]), row[1], row[2], [int(row[3]), int(row[4]), int(row[5])], toInt(row[7]))) for row in reader];
            
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
        #print i, levels[i]
        while slf.levels[i] > level:
            i = slf.parents[i];
            #print i, levels[i]
        return i;
            

Label = LabelInfo();
Label.initialize();



def labelPoints(points, labelimage, level = None):
    
    #points are (y,x,z) -> which is also the way the labeled image is read in
    x = points[:,0];
    y = points[:,1];
    z = points[:,2];    
    nPoint = x.size;    
    
    pointLabels = numpy.ones(nPoint, 'int32');
    if isinstance(labelimage, basestring):
        labelImage = io.readData(labelimage);
                  
    dsize = labelImage.shape;
    for i in range(nPoint):
        if y[i] >= 0 and y[i] < dsize[0] and x[i] >= 0 and x[i] < dsize[1] and z[i] >= 0 and z[i] < dsize[2]:
             pointLabels[i] = labelImage[y[i], x[i], z[i]];
    
    pointLabels = self.labelAtLevel(pointLabels, level);
    
    return pointLabels;

 
def labelAtLevel(label, level):
    global Label;
    
    if level is None:
        return label;
    else:
        return [Label.toLabelAtLevel(x, level) for x in label];

 
def countPointsInRegions(points, labelimage, level= None):
    pointLabels = self.labelPoints(points, labelimage, level = level); 
    cc = numpy.bincount(pointLabels);
    ii = numpy.nonzero(cc)[0];
    return numpy.vstack((ii,cc[ii])).T;
    
    
def labelToName(label):
    global Label;
    return [Label.name(x) for x in label];
    
def labelToAcronym(label):
    global Label;
    return [Label.acronym(x) for x in label];

def labelToColor(label):
    global Label;
    return [Label.color(x) for x in label];
 




def writePal(filename, cols):
    with open(filename, 'w') as f:
        for c in cols:
            f.write(('%3.3f' % c[0]).rjust(10) +  ('%3.3f' % c[1]).rjust(11) + ('%3.3f' % c[2]).rjust(11));
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
        writePal(filename, colarray);
        return colarray;
            
   
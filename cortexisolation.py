# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:17:38 2015

@author: mtllab
"""



import iDISCO.Analysis.Statistics as stat
import iDISCO.Analysis.Label as lbl
import iDISCO.IO.IO as io
import numpy



annotationFile = '/home/mtllab/Documents/warping/annotation_25_right_reslice.tif';
label = io.readData(annotationFile);
label = label.astype('int32');

labelids = numpy.unique(label);

outside = label == 0;
#for l in labelids:
#    if not lbl.labelAtLevel(l, 3) == 688: # remove everything except cortex
#        outside = numpy.logical_or(outside, label == l);
#for l in labelids:
#    if lbl.labelAtLevel(l, 5) == 1089: # remove hippocampus
#        outside = numpy.logical_or(outside, label == l);
#        
for l in labelids:
    if not lbl.labelAtLevel(l, 5) == 315: # remove everything except cortex
        outside = numpy.logical_or(outside, label == l);




############ These are to remove the structures curved around the lateral cortical plate        
#        
#for l in labelids:
#    if lbl.labelAtLevel(l, 6) == 31: # remove cingulate cortex
#        outside = numpy.logical_or(outside, label == l);
#        
#for l in labelids:
#    if lbl.labelAtLevel(l, 6) == 254: # remove retrosplenial cortex
#        outside = numpy.logical_or(outside, label == l);
#        
#for l in labelids:
#    if lbl.labelAtLevel(l, 6) == 631: # remove cortical Amygdala
#        outside = numpy.logical_or(outside, label == l);
#        
#for l in labelids:
#    if lbl.labelAtLevel(l, 5) == 698: # remove olfactory areas
#        outside = numpy.logical_or(outside, label == l);        
#
#for l in labelids:
#    if lbl.labelAtLevel(l, 4) == 703: # remove cortical subplate - amygdala
#        outside = numpy.logical_or(outside, label == l);       
#
#for l in labelids:
#    if lbl.labelAtLevel(l, 6) == 972: # remove prelimbic cortex
#        outside = numpy.logical_or(outside, label == l);    
        
#for l in labelids:
#    if not lbl.labelAtLevel(l, 6) == 453 : #somatosensory areas
#        outside = numpy.logical_or(outside, label == l);        
#for l in labelids:
#    if not lbl.labelAtLevel(l, 6) == 95: # keep insula
#        outside = numpy.logical_or(outside, label == l);
#    
#for l in labelids:
#    if not lbl.labelAtLevel(l, 6) == 961: # keep Piriform 
#        outside = numpy.logical_or(outside, label == l);
#
#for l in labelids:
#    if not lbl.labelAtLevel(l, 6) == 247: # keep auditory
#        outside = numpy.logical_or(outside, label == l);


cortex = label;
cortex[outside] = 0

io.writeData('/home/mtllab/Documents/cortex.mhd', cortex);



      
heatmap = io.readData('/home/mtllab/Documents/whiskers/exploration/mean_intact_allcells.mhd');

heatmap[outside] = 0;

io.writeData('/home/mtllab/Documents/whiskers/exploration/mean_intact_allcells_auditory.mhd', heatmap);

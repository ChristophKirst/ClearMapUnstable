# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 23:16:28 2015

@author: mtllab
"""



import iDISCO.Analysis.Statistics as stat
import iDISCO.Analysis.Label as lbl
import iDISCO.IO.IO as io
import numpy

group1 = ['/home/mtllab/Documents/whiskers/nonshaved/150620-6R/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/nonshaved/150620-7R/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/nonshaved/150620-9R/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/nonshaved/150620-10R/cells_heatmap.tif'];
                  
          
group2 = ['/home/mtllab/Documents/whisker806s/shaved/150620-1R/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/shaved/150620-2R/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/shaved/150620-3R/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/shaved/150620-4R/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/shaved/150620-5R/cells_heatmap.tif'];

g1 = stat.readGroup(group1);
g2 = stat.readGroup(group2);
g1.shape
g2.shape
806

pvals = stat.tTest(g1.astype('float'), g2.astype('float'));

pi = numpy.isnan(pvals);
pvals2 = pvals;
pvals2[pi] = 1.0;


pcutoff = 0.01;

pi = pvals2 > pcutoff;
pvals2[pi] = pcutoff;

annotationFile = '/home/mtllab/Documents/warping/annotation_25_right.tif';
label = io.readData(annotationFile);
label = label.astype('int32');

labelids = numpy.unique(label);

outside = label == 0;
for l in labelids:
    if lbl.labelAtLevel(l, 1) == 73: # ventricular system
        outside = numpy.logical_or(outside, label == l);

#io.writeData('/home/mtllab/Documents/whiskers/outisde.raw', outside.astype('int16'));

pvals2[outside] = pcutoff;

io.writeData('/home/mtllab/Documents/whiskers/pvalshaved.raw', pvals2);





# average and variance

g1 = stat.readGroup(group1);
g2 = stat.readGroup(group2);


g1a = numpy.mean(g1,axis = 0);
g1s = numpy.std(g1,axis = 0);

g2a = numpy.mean(g2,axis = 0);
g2s = numpy.std(g2,axis = 0);

io.writeData('/home/mtllab/Documents/whiskers/nonshaved_mean.raw', g1a);
io.writeData('/home/mtllab/Documents/whiskers/nonshaved_std.raw', g1s);

io.writeData('/home/mtllab/Documents/whiskers/shaved_mean.raw', g2a);
io.writeData('/home/mtllab/Documents/whiskers/shaved_std.raw', g2s);


#
#g1f = g1.astype('float');
#g2f = g2.astype('float');
#
# stat.tTest(g1f[:,:,210,80],  g2f[:,:,210,80])
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 23:16:28 2015

@author: mtllab
"""



import iDISCO.Analysis.Statistics as stat
import iDISCO.Analysis.Label as lbl
import iDISCO.IO.IO as io
import numpy, os

baseDirectory = '/home/mtllab/Documents/parental'


###########################################################################
################ Heat maps ################################################

group1 = ['/home/mtllab/Documents/parental/1Galcre/cells_heatmap.tif',
          '/home/mtllab/Documents/parental/5Galcre/cells_heatmap.tif'];
          
                  
group2 = ['/home/mtllab/Documents/parental/id2Ucn3cre/cells_heatmap.tif',
          '/home/mtllab/Documents/parental/id4Ucn3cre/cells_heatmap.tif'];

g1 = stat.readDataGroup(group1);
g2 = stat.readDataGroup(group2);
g1.shape
g2.shape


pvals, psign = stat.tTestVoxelization(g1.astype('float'), g2.astype('float'), signed = True, pcutoff = 0.05);
#pvals2 = stat.cutoffPValues(pvals, pcutoff = 0.05);


pvalsc = stat.colorPValues(pvals, psign, positive = [1,0], negative = [0,1]);


annotationFile = '/home/mtllab/Documents/warping/annotation_25_right_fullWD.tif';
label = io.readData(annotationFile);
label = label.astype('int32');

labelids = numpy.unique(label);

outside = label == 0;

for l in labelids:
    if lbl.labelAtLevel(l, 1) == 73: # remove ventricles
        outside = numpy.logical_or(outside, label == l);

pvalsc[outside]=0;

io.writeData(os.path.join(baseDirectory, 'pvalues_allcells.tif'), io.sagittalToCoronalData(pvalsc.astype('float32')));


g1a = numpy.mean(g1,axis = 0);
g1s = numpy.std(g1,axis = 0);

g2a = numpy.mean(g2,axis = 0);
g2s = numpy.std(g2,axis = 0);

io.writeData(os.path.join(baseDirectory, 'mean_females.raw'), io.sagittalToCoronalData(g1a));
io.writeData(os.path.join(baseDirectory, 'std_males.raw'), io.sagittalToCoronalData(g1s));

io.writeData(os.path.join(baseDirectory, 'mean_females.raw'), io.sagittalToCoronalData(g2a));
io.writeData(os.path.join(baseDirectory, 'std_males.raw'), io.sagittalToCoronalData(g2s));

##############################################################################
########## weighted Heat maps ################################################


group1 = ['/home/mtllab/Documents/parental/1Galcre/cells_heatmap_weighted.tif',
          '/home/mtllab/Documents/parental/5Galcre/cells_heatmap_weighted.tif'];
          
                  
group2 = ['/home/mtllab/Documents/parental/id2Ucn3cre/cells_heatmap_weighted.tif',
          '/home/mtllab/Documents/parental/id4Ucn3cre/cells_heatmap_weighted.tif'];


g1 = stat.readDataGroup(group1);
g2 = stat.readDataGroup(group2);
g1.shape
g2.shape


pvals, psign = stat.tTestVoxelization(g1.astype('float'), g2.astype('float'), signed = True, pcutoff = 0.05);
#pvals2 = stat.cutoffPValues(pvals, pcutoff = 0.05);


pvalsc = stat.colorPValues(pvals, psign, positive = [1,0], negative = [0,1]);


annotationFile = '/home/mtllab/Documents/warping/annotation_25_right_fullWD.tif';
label = io.readData(annotationFile);
label = label.astype('int32');

labelids = numpy.unique(label);

outside = label == 0;

for l in labelids:
    if lbl.labelAtLevel(l, 1) == 73: # remove ventricles
        outside = numpy.logical_or(outside, label == l);

pvalsc[outside]=0;

io.writeData(os.path.join(baseDirectory, 'pvalues_weighted.tif'), io.sagittalToCoronalData(pvalsc.astype('float32')));

g1a = numpy.mean(g1,axis = 0);
g1s = numpy.std(g1,axis = 0);

g2a = numpy.mean(g2,axis = 0);
g2s = numpy.std(g2,axis = 0);

io.writeData(os.path.join(baseDirectory, 'mean_females_weighted.raw'), io.sagittalToCoronalData(g1a));
io.writeData(os.path.join(baseDirectory, 'std_females_weighted.raw'), io.sagittalToCoronalData(g1s));

io.writeData(os.path.join(baseDirectory, 'mean_males_weighted.raw'), io.sagittalToCoronalData(g2a));
io.writeData(os.path.join(baseDirectory, 'std_males_weighted.raw'), io.sagittalToCoronalData(g2s));





###################################################################################
################# Cell counting per regions #######################################


group1 = ['/home/mtllab/Documents/parental/1Galcre/cells_transformed_to_Atlas.npy',
          '/home/mtllab/Documents/parental/5Galcre/cells_transformed_to_Atlas.npy'];
          
                  
group2 = ['/home/mtllab/Documents/parental/id2Ucn3cre/cells_transformed_to_Atlas.npy',
          '/home/mtllab/Documents/parental/id4Ucn3cre/cells_transformed_to_Atlas.npy'];
 


ids, pc1 = stat.countPointsGroupInRegions(group1, withIds = True, labeledImage = annotationFile);
pc2 = stat.countPointsGroupInRegions(group2, withIds = False, labeledImage = annotationFile);

pvals, psign = stat.tTestPointsInRegions(pc1, pc2, pcutoff = None, signed = True);


annotationFile = '/home/mtllab/Documents/warping/annotation_25_right_fullWD.tif';
label = io.readData(annotationFile);
label = label.astype('int32');

labelids = numpy.unique(label);

outside = label == 0;

for l in labelids:
    if lbl.labelAtLevel(l, 1) == 73: # remove ventricles
        outside = numpy.logical_or(outside, label == l);

pvals[outside]=0;
#make table


table = numpy.zeros(ids.shape, dtype=[('id','int64'),('mean1','f8'),('mean2','f8'),('pvalue', 'f8'),('psign', 'int64'),('name', 'a256')])
table["id"] = ids;
table["mean1"] = pc1.mean(axis = 1);
table["mean2"] = pc2.mean(axis = 1);
table["pvalue"] = pvals;
table["psign"] = psign;
table["name"] = lbl.labelToName(ids);


#sort by pvalue
ii = numpy.argsort(pvals);
tableSorted = table.copy();
tableSorted = tableSorted[ii];

with open(os.path.join(baseDirectory, 'regions_pvalues.csv'),'w') as f:
    for sublist in tableSorted:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();



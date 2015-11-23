# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 23:16:28 2015

@author: mtllab
"""



import iDISCO.Analysis.Statistics as stat
import iDISCO.Analysis.Label as lbl
import iDISCO.IO.IO as io
import numpy, os

baseDirectory = '/home/mtllab/Documents/whiskers/exploration/'


###########################################################################
################ Heat maps ################################################

group1 = [#'/home/mtllab/Documents/whiskers/exploration/1/cells_heatmap.tif',
          #'/home/mtllab/Documents/whiskers/exploration/2/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/exploration/3/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/exploration/4/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/exploration/5/cells_heatmap.tif'];
          
                  
group2 = ['/home/mtllab/Documents/whiskers/exploration/6/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/exploration/7/cells_heatmap.tif',
         # '/home/mtllab/Documents/whiskers/exploration/8/oldcounts/cells_heatmap.tif',
          '/home/mtllab/Documents/whiskers/exploration/9/cells_heatmap.tif']
         # '/home/mtllab/Documents/whiskers/exploration/10/cells_heatmap.tif'];

g1 = stat.readDataGroup(group1);
g2 = stat.readDataGroup(group2);
g1.shape
g2.shape


pvals, psign = stat.tTestVoxelization(g1.astype('float'), g2.astype('float'), signed = True, pcutoff = 0.05);
#pvals2 = stat.cutoffPValues(pvals, pcutoff = 0.05);


pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);


annotationFile = '/home/mtllab/Documents/warping/annotation_25_right.tif';
label = io.readData(annotationFile);
label = label.astype('int32');

labelids = numpy.unique(label);

outside = label == 0;
#
#for l in labelids:
#    if lbl.labelAtLevel(l, 1) == 73: # remove ventricles
#        outside = numpy.logical_or(outside, label == l);
#
pvalsc[outside]=0;

io.writeData(os.path.join(baseDirectory, 'pvalues_allcells.tif'), io.sagittalToCoronalData(pvalsc.astype('float32')));


g1a = numpy.mean(g1,axis = 0);
g1s = numpy.std(g1,axis = 0);

g2a = numpy.mean(g2,axis = 0);
g2s = numpy.std(g2,axis = 0);
#
g1a[outside] = 0;
g1s[outside] = 0;
g2a[outside] = 0;
g2s[outside] = 0;


io.writeData(os.path.join(baseDirectory, 'mean_shaved_allcells.raw'), io.sagittalToCoronalData(g1a));
io.writeData(os.path.join(baseDirectory, 'std_shaved_allcells.raw'), io.sagittalToCoronalData(g1s));

io.writeData(os.path.join(baseDirectory, 'mean_intact_allcells.raw'), io.sagittalToCoronalData(g2a));
io.writeData(os.path.join(baseDirectory, 'std_intact_allcells.raw'), io.sagittalToCoronalData(g2s));

##############################################################################
########## weighted Heat maps ################################################


import iDISCO.Analysis.Statistics as stat
import iDISCO.Analysis.Label as lbl
import iDISCO.IO.IO as io
import numpy, os

baseDirectory = '/home/mtllab/Documents/whiskers/exploration/'



group1 = [#'/home/mtllab/Documents/whiskers/exploration/1/cells_heatmap_weighted.tif',
          #'/home/mtllab/Documents/whiskers/exploration/2/cells_heatmap_weighted.tif',
          '/home/mtllab/Documents/whiskers/exploration/3/cells_heatmap_weighted.tif',
          '/home/mtllab/Documents/whiskers/exploration/4/cells_heatmap_weighted.tif',
          '/home/mtllab/Documents/whiskers/exploration/5/cells_heatmap_weighted.tif'];
          
                  
group2 = ['/home/mtllab/Documents/whiskers/exploration/6/cells_heatmap_weighted.tif',
          '/home/mtllab/Documents/whiskers/exploration/7/cells_heatmap_weighted.tif',
       #   '/home/mtllab/Documents/whiskers/exploration/8/oldcounts/cells_heatmap_weighted.tif',
          '/home/mtllab/Documents/whiskers/exploration/9/cells_heatmap_weighted.tif']
       #   '/home/mtllab/Documents/whiskers/exploration/10/oldcounts/cells_heatmap_weighted.tif'];

g1 = stat.readDataGroup(group1);
g2 = stat.readDataGroup(group2);
g1.shape
g2.shape


pvals, psign = stat.tTestVoxelization(g1.astype('float'), g2.astype('float'), signed = True, pcutoff = 0.05);
#pvals2 = stat.cutoffPValues(pvals, pcutoff = 0.05);

annotationFile = '/home/mtllab/Documents/warping/annotation_25_right.tif';
label = io.readData(annotationFile);
label = label.astype('int32');

labelids = numpy.unique(label);

outside = label == 0;

for l in labelids:
    if lbl.labelAtLevel(l, 1) == 73: # remove ventricles
        outside = numpy.logical_or(outside, label == l);

pvalsc[outside]=0;


pvalsc = stat.colorPValues(pvals, psign, positive = [0,1], negative = [1,0]);

io.writeData(os.path.join(baseDirectory, 'pvalues_weighted.tif'), io.sagittalToCoronalData(pvalsc.astype('float32')));

g1a = numpy.mean(g1,axis = 0);
g1s = numpy.std(g1,axis = 0);

g2a = numpy.mean(g2,axis = 0);
g2s = numpy.std(g2,axis = 0);

io.writeData(os.path.join(baseDirectory, 'mean_shaved_weighted.raw'), io.sagittalToCoronalData(g1a));
io.writeData(os.path.join(baseDirectory, 'std_shaved_weighted.raw'), io.sagittalToCoronalData(g1s));

io.writeData(os.path.join(baseDirectory, 'mean_intact_weighted.raw'), io.sagittalToCoronalData(g2a));
io.writeData(os.path.join(baseDirectory, 'std_intact_weighted.raw'), io.sagittalToCoronalData(g2s));




###################################################################################
################# Cell counting per regions #######################################


import iDISCO.Analysis.Statistics as stat
import iDISCO.Analysis.Label as lbl
import iDISCO.IO.IO as io
import numpy, os

baseDirectory = '/home/mtllab/Documents/whiskers/exploration/'



group1 = [#'/home/mtllab/Documents/whiskers/exploration/1/cells_transformed_to_Atlas.npy',
          #'/home/mtllab/Documents/whiskers/exploration/2/cells_transformed_to_Atlas.npy',
          '/home/mtllab/Documents/whiskers/exploration/3/cells_transformed_to_Atlas.npy',
          '/home/mtllab/Documents/whiskers/exploration/4/cells_transformed_to_Atlas.npy',
          '/home/mtllab/Documents/whiskers/exploration/5/cells_transformed_to_Atlas.npy'];
          
                  
group2 = ['/home/mtllab/Documents/whiskers/exploration/6/cells_transformed_to_Atlas.npy',
          '/home/mtllab/Documents/whiskers/exploration/7/cells_transformed_to_Atlas.npy',
         # '/home/mtllab/Documents/whiskers/exploration/8/cells_transformed_to_Atlas.npy',
          '/home/mtllab/Documents/whiskers/exploration/9/cells_transformed_to_Atlas.npy']
         # '/home/mtllab/Documents/whiskers/exploration/10/cells_transformed_to_Atlas.npy'];

group1i = [fn.replace('cells_transformed_to_Atlas', 'intensities') for fn in group1];

group2i = [fn.replace('cells_transformed_to_Atlas', 'intensities') for fn in group2];


ABAlabeledImage = '/home/mtllab/Documents/warping/annotation_25_right.tif';


ids, pc1, pc1i = stat.countPointsGroupInRegions(group1, intensityGroup = group1i, returnIds = True, labeledImage = ABAlabeledImage, returnCounts = True, collapse = 'yes');
pc2, pc2i = stat.countPointsGroupInRegions(group2, intensityGroup = group2i, returnIds = False, labeledImage = ABAlabeledImage, returnCounts = True, collapse = 'yes');


pvals, psign = stat.tTestPointsInRegions(pc1, pc2, pcutoff = None, signed = True, equal_var = True);
pvalsi, psigni = stat.tTestPointsInRegions(pc1i, pc2i, pcutoff = None, signed = True, equal_var = True);


import iDISCO.Analysis.Tools.MultipleComparisonCorrection as mc

lowcount = numpy.sum(pc1, axis=1) + numpy.sum(pc2, axis=1)
iid = lowcount > 200;

ids0 = ids[iid];
pc1i0 = pc1i[iid];
pc2i0 = pc2i[iid];
pc10 = pc1[iid];
pc20 = pc2[iid];
psigni0 = psigni[iid];
psign0 = psign[iid];
pvalsi0 = pvalsi[iid];
pvals0 = pvals[iid];
qvalsi0 = mc.estimateQValues(pvalsi0);
qvals0 = mc.estimateQValues(pvals0);

pvalsiBH0 = mc.correctPValues(pvalsi0);
pvalsBH0 = mc.correctPValues(pvals0);


#make table

dtypes = [('id','int64'),('mean1','f8'),('std1','f8'),('mean2','f8'),('std2','f8'),('pvalue', 'f8'),('pvalue-BH', 'f8'),('qvalue', 'f8'),('psign', 'int64')];
for i in range(len(group1)):
    dtypes.append(('count1_%d' % i, 'f8'));
for i in range(len(group2)):
    dtypes.append(('count2_%d' % i, 'f8'));   
dtypes.append(('name', 'a256'));

table = numpy.zeros(ids0.shape, dtype = dtypes)
table["id"] = ids0;
table["mean1"] = pc1i0.mean(axis = 1)/1000;
table["std1"] = pc1i0.std(axis = 1)/1000;
table["mean2"] = pc2i0.mean(axis = 1)/1000;
table["std2"] = pc2i0.std(axis = 1)/1000;
table["pvalue"] = pvalsi0;
table["pvalue-BH"] = pvalsiBH0;
table["qvalue"] = qvalsi0;

table["psign"] = psigni0;
for i in range(len(group1)):
    table["count1_%d" % i] = pc10[:,i];
for i in range(len(group2)):
    table["count2_%d" % i] = pc20[:,i];
table["name"] = lbl.labelToName(ids0);


#sort by pvalue
ii = numpy.argsort(pvalsi0);
tableSorted = table.copy();
tableSorted = tableSorted[ii];

with open(os.path.join(baseDirectory, 'qvalues_intensities_n4_collapsedregions_test.csv'),'w') as f:
    f.write(', '.join([str(item) for item in table.dtype.names]));
    f.write('\n');
    for sublist in tableSorted:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

#############################


#make table

dtypes = [('id','int64'),('mean1','f8'),('std1','f8'),('mean2','f8'),('std2','f8'),('pvalue', 'f8'),('pvalue-BH', 'f8'),('qvalue', 'f8'),('psign', 'int64')];
for i in range(len(group1)):
    dtypes.append(('count1_%d' % i, 'f8'));
for i in range(len(group2)):
    dtypes.append(('count2_%d' % i, 'f8'));   
dtypes.append(('name', 'a256'));

table = numpy.zeros(ids0.shape, dtype = dtypes)
table["id"] = ids0;
table["mean1"] = pc10.mean(axis = 1);
table["std1"] = pc10.std(axis = 1);
table["mean2"] = pc20.mean(axis = 1);
table["std2"] = pc20.std(axis = 1);
table["pvalue"] = pvals0;
table["pvalue-BH"] = pvalsBH0;
table["qvalue"] = qvals0;
table["psign"] = psign0;
for i in range(len(group1)):
    table["count1_%d" % i] = pc10[:,i];
for i in range(len(group2)):
    table["count2_%d" % i] = pc20[:,i];
table["name"] = lbl.labelToName(ids0);


#sort by pvalue
ii = numpy.argsort(pvals0);
tableSorted = table.copy();
tableSorted = tableSorted[ii];

with open(os.path.join(baseDirectory, 'qvalues_counts_n4_collapsedregions_test.csv'),'w') as f:
    f.write(', '.join([str(item) for item in table.dtype.names]));
    f.write('\n');
    for sublist in tableSorted:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();



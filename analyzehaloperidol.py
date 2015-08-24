# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 23:16:28 2015

@author: mtllab
"""



import iDISCO.Analysis.Statistics as stat
import iDISCO.Analysis.Label as lbl
import iDISCO.IO.IO as io
import numpy, os

baseDirectory = '/home/mtllab/Documents/Haloperidol'


group1 = ['/home/mtllab/Documents/Haloperidol/1266/cells_heatmap.tif',
          '/home/mtllab/Documents/Haloperidol/1267/cells_heatmap.tif',
          '/home/mtllab/Documents/Haloperidol/1268/cells_heatmap.tif',
          '/home/mtllab/Documents/Haloperidol/1269/cells_heatmap.tif',
          '/home/mtllab/Documents/Haloperidol/1270/cells_heatmap.tif'];
          
          
          
group2 = ['/home/mtllab/Documents/Haloperidol/1271/cells_heatmap.tif',
          '/home/mtllab/Documents/Haloperidol/1272/cells_heatmap.tif',
          '/home/mtllab/Documents/Haloperidol/1273/cells_heatmap.tif',
          '/home/mtllab/Documents/Haloperidol/1274/cells_heatmap.tif',
          '/home/mtllab/Documents/Haloperidol/1275/cells_heatmap.tif'];

g1 = stat.readGroup(group1);
g2 = stat.readGroup(group2);
g1.shape
g2.shape


pvals = stat.tTest(g1.astype('float'), g2.astype('float'), signed = False);

pi = numpy.isnan(pvals);
pvals2 = pvals.copy();
pvals2[pi] = 1.0;

pcutoff = 0.05;

pi = pvals2 > pcutoff;
pvals2[pi] = pcutoff;
pi = pvals2 < - pcutoff;
pvals2[pi] = - pcutoff;


io.writeData(os.path.join(baseDirectory, 'pvalues_allcells.raw'), io.sagittalToCoronalData(numpy.abs(pvals2)));


pi = pvals2 > pcutoff;
pvals2[pi] = pcutoff;
pi = pvals2 < - pcutoff;
pvals2[pi] = - pcutoff;




pvalsColor = stat.colorPValues(1 / pvals2 * pcutoff, positive = [1,0,0], negative = [0,1,0]);
pvalsColor = pvalsColor[:,:,:,0:2];
pvalsColor = 255*255 * pvalsColor;
#pvalsColor = pvalsColor.astype('int32');
pvalsColor = pvalsColor.astype('float32');
pvalsColor.max()
pvalsColor.min();

io.writeData(os.path.join(baseDirectory, 'pvalues.tif'), pvalsColor);


#annotationFile = '/home/mtllab/Documents/warping/annotation_25_right.tif';
#label = io.readData(annotationFile);
#label = label.astype('int32');
#
#labelids = numpy.unique(label);
#
#outside = label == 0;
#for l in labelids:
#    if lbl.labelAtLevel(l, 1) == 73: # ventricular system
#        outside = numpy.logical_or(outside, label == l);
#
##io.writeData('/home/mtllab/Documents/outisde.raw', outside.astype('int16'));
#
#pvals2[outside] = pcutoff;
#
#io.writeData(os.path.join(baseDirectory, 'pval_haloperidol_vs_saline_p005.raw'), pvals2);




# average and variance

g1 = stat.readGroup(group1);
g2 = stat.readGroup(group2);


g1a = numpy.mean(g1,axis = 0);
g1s = numpy.std(g1,axis = 0);

g2a = numpy.mean(g2,axis = 0);
g2s = numpy.std(g2,axis = 0);

io.writeData(os.path.join(baseDirectory, 'mean_haloperidol.raw'), g1a);
io.writeData(os.path.join(baseDirectory, 'std_haloperidol.raw'), g1s);

io.writeData(os.path.join(baseDirectory, 'mean_saline.raw'), g2a);
io.writeData(os.path.join(baseDirectory, 'std_saline.raw'), g2s);


#
#g1f = g1.astype('float');
#g2f = g2.astype('float');
#
# stat.tTest(g1f[:,:,210,80],  g2f[:,:,210,80])





##################################################
io.sagittalToCoronalData('/home/mtllab/Documents/Haloperidol/1271/cells_heatmap.tif', '/home/mtllab/Documents/Haloperidol/test.tif')



from iDISCO.Analysis.Voxelization import voxelize
import iDISCO.IO.IO as io

points = io.readPoints('/home/mtllab/Documents/Haloperidol/1270/cells_transformed_to_Atlas.csv');
points = points[:,[1,0,2]];

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1270/intensities.csv');

#intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1270/cells_heatmap.tif', voxelizationSize = (15,15,15), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1270/cells_intensity_heatmap.tif', io.sagittalToCoronalData(voximg.astype('float32')))

#############################


from iDISCO.Analysis.Voxelization import voxelize
import iDISCO.IO.IO as io

points = io.readPoints('/home/mtllab/Documents/Haloperidol/1270/cells_transformed_to_Atlas.csv');
points = points[:,[1,0,2]];

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1270/intensities.csv');

intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1268/cells_heatmap.tif', voxelizationSize = (5,5,5), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1270/cells_intensity_heatmap5.tif', io.sagittalToCoronalData(voximg.astype('float32')))

##################################


from iDISCO.Analysis.Voxelization import voxelize
import iDISCO.IO.IO as io
import numpy, os

directories = ['/home/mtllab/Documents/Haloperidol/1266/',
               '/home/mtllab/Documents/Haloperidol/1268/',
               '/home/mtllab/Documents/Haloperidol/1269/',
               '/home/mtllab/Documents/Haloperidol/1270/'];


pointFiles = [os.path.join(x,'cells_transformed_to_Atlas.csv') for x in directories];
intensityFiles = [os.path.join(x,'intensities.csv') for x in directories];

                  
points = [];
intensities = [];
for i in range(len(pointFiles)):
    pts = io.readPoints(pointFiles[i]);
    pts = pts[:,[1,0,2]];
    points.append(pts);
    
    its = io.readPoints(intensityFiles[i]);
    its = its[its > 0];
    intensities.append(its);

points = numpy.concatenate(points);
intensities = numpy.concatenate(points);

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1270/cells_heatmap.tif', voxelizationSize = (5,5,5), voxelizationWeights = intensities);
io.writeData('/home/mtllab/Documents/Haloperidol/1267/cells_intensity_heatmap5.tif', io.sagittalToCoronalData(voximg.astype('float32')))













points = io.readPoints('/home/mtllab/Documents/Haloperidol/1267/cells_transformed_to_Atlas.csv');



intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1267/intensities.csv');

#intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1270/cells_heatmap.tif', voxelizationSize = (5,5,5), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1267/cells_intensity_heatmap5.tif', io.sagittalToCoronalData(voximg.astype('float32')))

##################################





points = io.readPoints('/home/mtllab/Documents/Haloperidol/1271/cells.csv');
intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1271/intensities.csv');


voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1271/cells_heatmap.tif', voxelizationSize = (15,15,15), voxelizationWeights = intensities);


#####################################


import iDISCO.Analysis.Statistics as stat
import iDISCO.Analysis.Label as lbl
import iDISCO.IO.IO as io
import numpy, os

baseDirectory = '/home/mtllab/Documents/Haloperidol'


group1 = ['/home/mtllab/Documents/Haloperidol/1266/cells_intensity_heatmap10.tif',
          '/home/mtllab/Documents/Haloperidol/1267/cells_intensity_heatmap10.tif',
          '/home/mtllab/Documents/Haloperidol/1268/cells_intensity_heatmap10.tif',
          '/home/mtllab/Documents/Haloperidol/1269/cells_intensity_heatmap10.tif',
          '/home/mtllab/Documents/Haloperidol/1270/cells_intensity_heatmap10.tif'];
          
          
          
group2 = ['/home/mtllab/Documents/Haloperidol/1271/cells_intensity_heatmap10.tif',
          '/home/mtllab/Documents/Haloperidol/1272/cells_intensity_heatmap10.tif',
          '/home/mtllab/Documents/Haloperidol/1273/cells_intensity_heatmap10.tif',
          '/home/mtllab/Documents/Haloperidol/1274/cells_intensity_heatmap10.tif',
          '/home/mtllab/Documents/Haloperidol/1275/cells_intensity_heatmap10.tif'];

g1 = stat.readGroup(group1);
g2 = stat.readGroup(group2);
g1.shape
g2.shape


pvals = stat.tTest(g1.astype('float'), g2.astype('float'), signed = False);

pi = numpy.isnan(pvals);
pvals2 = pvals.copy();
pvals2[pi] = 1.0;

pcutoff = 0.05;

pi = pvals2 > pcutoff;
pvals2[pi] = pcutoff;
#pi = pvals2 < - pcutoff;
#pvals2[pi] = - pcutoff;


io.writeData(os.path.join(baseDirectory, 'pvalues_intensitiesweighted10.raw'), numpy.abs(pvals2));



#pvalsColor = stat.colorPValues(1 / pvals2 * pcutoff, positive = [1,0,0], negative = [0,1,0]);
#pvalsColor = pvalsColor[:,:,:,0:2];
#pvalsColor = 255*255 * pvalsColor;
#pvalsColor = pvalsColor.astype('int32');
#pvalsColor = pvalsColor.astype('float32');
#pvalsColor.max()
#pvalsColor.min();

#io.writeData(os.path.join(baseDirectory, 'pvalues.tif'), pvalsColor);


#annotationFile = '/home/mtllab/Documents/warping/annotation_25_right.tif';
#label = io.readData(annotationFile);
#label = label.astype('int32');
#
#labelids = numpy.unique(label);
#
#outside = label == 0;
#for l in labelids:
#    if lbl.labelAtLevel(l, 1) == 73: # ventricular system
#        outside = numpy.logical_or(outside, label == l);
#
##io.writeData('/home/mtllab/Documents/outisde.raw', outside.astype('int16'));
#
#pvals2[outside] = pcutoff;
#
#io.writeData(os.path.join(baseDirectory, 'pval_haloperidol_vs_saline_p005.raw'), pvals2);




# average and variance

g1 = stat.readGroup(group1);
g2 = stat.readGroup(group2);


g1a = numpy.mean(g1,axis = 0);
g1s = numpy.std(g1,axis = 0);

g2a = numpy.mean(g2,axis = 0);
g2s = numpy.std(g2,axis = 0);

io.writeData(os.path.join(baseDirectory, 'mean_haloperidol_int10.raw'), g1a);
io.writeData(os.path.join(baseDirectory, 'std_haloperidol_int10.raw'), g1s);

io.writeData(os.path.join(baseDirectory, 'mean_saline_int10.raw'), g2a);
io.writeData(os.path.join(baseDirectory, 'std_saline_int10.raw'), g2s);


#
#g1f = g1.astype('float');
#g2f = g2.astype('float');
#
# stat.tTest(g1f[:,:,210,80],  g2f[:,:,210,80])

###################################################################

from iDISCO.Analysis.Voxelization import voxelize
import iDISCO.IO.IO as io

points = io.readPoints('/home/mtllab/Documents/Haloperidol/1266/cells_transformed_to_Atlas.csv');
points = points[:,[1,0,2]];

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1266/intensities.csv');

intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1268/cells_heatmap.tif', voxelizationSize = (10,10,10), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1266/cells_intensity_heatmap10.tif', io.sagittalToCoronalData(voximg.astype('float32')))









points = io.readPoints('/home/mtllab/Documents/Haloperidol/1267/cells_transformed_to_Atlas.csv');

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1267/intensities.csv');

intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1268/cells_heatmap.tif', voxelizationSize = (10,10,10), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1267/cells_intensity_heatmap10.tif', io.sagittalToCoronalData(voximg.astype('float32')))







points = io.readPoints('/home/mtllab/Documents/Haloperidol/1268/cells_transformed_to_Atlas.csv');
points = points[:,[1,0,2]];

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1268/intensities.csv');

intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1268/cells_heatmap.tif', voxelizationSize = (10,10,10), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1268/cells_intensity_heatmap10.tif', io.sagittalToCoronalData(voximg.astype('float32')))








points = io.readPoints('/home/mtllab/Documents/Haloperidol/1269/cells_transformed_to_Atlas.csv');
points = points[:,[1,0,2]];

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1269/intensities.csv');

intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1269/cells_heatmap.tif', voxelizationSize = (10,10,10), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1269/cells_intensity_heatmap10.tif', io.sagittalToCoronalData(voximg.astype('float32')))







points = io.readPoints('/home/mtllab/Documents/Haloperidol/1270/cells_transformed_to_Atlas.csv');
points = points[:,[1,0,2]];

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1266/intensities.csv');

intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1270/cells_heatmap.tif', voxelizationSize = (10,10,10), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1270/cells_intensity_heatmap10.tif', io.sagittalToCoronalData(voximg.astype('float32')))








points = io.readPoints('/home/mtllab/Documents/Haloperidol/1271/cells_transformed_to_Atlas.csv');
points = points[:,[1,0,2]];

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1271/intensities.csv');

intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1268/cells_heatmap.tif', voxelizationSize = (10,10,10), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1271/cells_intensity_heatmap10.tif', io.sagittalToCoronalData(voximg.astype('float32')))






points = io.readPoints('/home/mtllab/Documents/Haloperidol/1272/cells_transformed_to_Atlas.csv');
points = points[:,[1,0,2]];

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1272/intensities.csv');

intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1268/cells_heatmap.tif', voxelizationSize = (10,10,10), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1272/cells_intensity_heatmap10.tif', io.sagittalToCoronalData(voximg.astype('float32')))





points = io.readPoints('/home/mtllab/Documents/Haloperidol/1273/cells_transformed_to_Atlas.csv');
points = points[:,[1,0,2]];

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1273/intensities.csv');

intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1268/cells_heatmap.tif', voxelizationSize = (10,10,10), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1273/cells_intensity_heatmap10.tif', io.sagittalToCoronalData(voximg.astype('float32')))




points = io.readPoints('/home/mtllab/Documents/Haloperidol/1274/cells_transformed_to_Atlas.csv');
points = points[:,[1,0,2]];

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1274/intensities.csv');

intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1268/cells_heatmap.tif', voxelizationSize = (10,10,10), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1274/cells_intensity_heatmap10.tif', io.sagittalToCoronalData(voximg.astype('float32')))



points = io.readPoints('/home/mtllab/Documents/Haloperidol/1275/cells_transformed_to_Atlas.csv');
points = points[:,[1,0,2]];

intensities = io.readPoints('/home/mtllab/Documents/Haloperidol/1275/intensities.csv');

intensities = intensities[intensities > 0];

points.shape
intensities.shape

voximg = voxelize(points, '/home/mtllab/Documents/Haloperidol/1275/cells_heatmap.tif', voxelizationSize = (10,10,10), voxelizationWeights = intensities);

io.writeData('/home/mtllab/Documents/Haloperidol/1275/cells_intensity_heatmap10.tif', io.sagittalToCoronalData(voximg.astype('float32')))

# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 16:03:57 2015

@author: mtllab
"""

import iDISCO.Visualization.Plot as plt;
import iDISCO.IO.IO as io;

import numpy

dataSource =  '/home/mtllab/Documents/Haloperidol/1268/150819_0_8X-cfos_14-22-33/14-22-33_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif'
pointSource= '/home/mtllab/Documents/Haloperidol/1268/cells.npy';
x = (100,1000);
y = (600,1500);
z = (300,1200);


data = plt.overlayPoints(dataSource, pointSource, x = x, y = y, z = z, pointColor = None);
io.writeData('/home/mtllab/Documents/Haloperidol/1268/cells_check.tif', data)


########################################
########################################
import iDISCO.IO.IO as io;
import numpy as np;
import matplotlib.pyplot as matpp;
import iDISCO.Analysis.Label as lbl;
import iDISCO.Analysis.Statistics as stats

annotationFile = '/home/mtllab/Documents/warping/annotation_25_right.tif';

fnsH = [('/home/mtllab/Documents/Haloperidol/%d' % xx) for xx in range(1266, 1271)];
fnsH = [('/home/mtllab/Documents/Haloperidol/%d' % xx) for xx in range(1266, 1267)];
pointsH = [];
intensitiesH = [];
for i in range(len(fnsH)):
    phalop, ihalop = io.readPoints(((fnsH[i] + '/cells_transformed_to_Atlas.npy'),  (fnsH[i] + '/intensities.npy')));
    pointsH.append(phalop);
    intensitiesH.append(ihalop); 


fnsS = [('/home/mtllab/Documents/Haloperidol/%d' % xx) for xx in range(1271, 1276)];
fnsS = [('/home/mtllab/Documents/Haloperidol/%d' % xx) for xx in range(1271, 1272)];
pointsS = [];
intensitiesS = [];
for i in range(len(fnsH)):
    psaline, isaline = io.readPoints(((fnsS[i] + '/cells_transformed_to_Atlas.npy'),  (fnsS[i] + '/intensities.npy')));
    pointsS.append(psaline);
    intensitiesS.append(isaline); 


phalop = np.concatenate(pointsH);
ihalop = np.concatenate(intensitiesH);

psaline = np.concatenate(pointsS);
isaline= np.concatenate(intensitiesS);




lblpointsH = lbl.labelPoints(phalop, annotationFile);
lblpointsS = lbl.labelPoints(psaline, annotationFile);



#(p,s,n1,n2) = stats.testCompletedCumulativesInSpheres(phalop, ihalop[:,3], psaline, isaline[:,3], dataSize = annotationFile, radius = 2);




# cumulative statistics in each region
label = io.readData(annotationFile);
label = label.astype('int32');

labelids = np.unique(label);

ii = 3;

pval = np.zeros(labelids.size);
sval = np.zeros(labelids.size);
dn = pval.copy();
di = pval.copy();

for i in range(labelids.size):
    
    print "%d / %d %s" % (i, len(labelids), lbl.labelToName([labelids[i]])[0]);
    l = labelids[i];
    
    iH = lblpointsH == l;
    iS = lblpointsS == l;
    
    nH = iH.sum();
    nS = iS.sum();
    
    dn[i] = nH - nS;
    
    inH = ihalop[iH, ii];
    inS = isaline[iS, ii];    
    di[i] = inH.sum() - inS.sum();

    if nH > 0 and nS > 0:    
    
        #pH = phalop[iH,:];
        #pS = psaline[iS,:];
 
        try:
            (pp, ss) = stats.testCompletedCumulatives((inH, inS), method = "KS");
        except:
            pp = 0;
            ss = 0;
        
        pval[i] = pp;
        sval[i] = ss;
    
    else:
        pval[i] = 1;
        sval[i] = 0;
    
        

pval

dtypes = [('id','int64'),('pvalue', 'f8'),('svalue', 'f8'),('dn', 'f8'),('di', 'f8'),('name', 'a256')];

ifilter = di > 5000;
fids = labelids[ifilter];
fpval = pval[ifilter];
fsval = sval[ifilter];
fdn = dn[ifilter];
fdi = di[ifilter];


table = np.zeros(fids.shape, dtype = dtypes)
table["id"] = fids;
table["pvalue"] = fpval;
table["svalue"] = fsval;
table["dn"] = fdn;
table["di"] = fdi;

table["name"] = lbl.labelToName(fids);


#sort by pvalue
ii = np.argsort(fpval);
tableSorted = table.copy();
tableSorted = tableSorted[ii];

import os;
with open(os.path.join('/home/mtllab/Documents/Haloperidol', 'pvalues-KS-regions.csv'),'w') as f:
    f.write(', '.join([str(item) for item in table.dtype.names]));
    f.write('\n');
    for sublist in tableSorted:
        f.write(', '.join([str(item) for item in sublist]));
        f.write('\n');
    f.close();

#outside = numpy.zeros(label.shape, dtype = bool);

#for l in labelids:
#    if not (lbl.labelAtLevel(l, 5) == 315): #or (lbl.labelAtLevel(l, 4) == 703) or (lbl.labelAtLevel(l, 5) == 698)):
#        outside = numpy.logical_or(outside, label == l);












l = 440;
ii = 3;
iH = lblpointsH == l;
iS = lblpointsS == l;

nH = iH.sum();
nS = iS.sum();

inH = ihalop[iH, ii];
inS = isaline[iS, ii];    
di = inH.sum() - inS.sum();

if nH > 0 and nS > 0: 
    plt.figure()    
    try:
        (pp, ss) = stats.testCompletedCumulatives((inH, inS), method = "AD", plot = True);
    except:
        pp = 0;
        ss = 0;

print pp, ss



lblpH = []; lblpS =[];
for i in range(len(pointsH)):
    lblpH.append(lbl.labelPoints(pointsH[i], annotationFile));
for i in range(len(pointsS)):
    lblpS.append(lbl.labelPoints(pointsS[i], annotationFile));



#plt.figure()

l = 1031;
l = 672;
ii = 3;
lHH = []; lSS = [];
iHH = []; iSS = [];
for i in range(len(lblpH)):
    lHH.append(lblpH[i] == l);
    print i
    if lHH[-1].sum() > 0:
        iHH.append(intensitiesH[i][lHH[-1],ii]);
for i in range(len(lblpS)):
    lSS.append(lblpS[i] == l);
    if lSS[-1].sum() > 0:
        iSS.append(intensitiesS[i][lSS[-1],ii]);

import matplotlib.pyplot as plt
plt.figure()
(pp, ss) = stats.testCompletedCumulatives(iHH, method = "AD", plot = True);

plt.figure()
(pp, ss) = stats.testCompletedCumulatives(iSS, method = "AD", plot = True);

    


(pp, ss) = stats.testCompletedInvertedCumulatives(iHH, method = "AD", plot = True);

#plt.figure()
(pp, ss) = stats.testCompletedInvertedCumulatives(iSS, method = "AD", plot = True);




#plt.figure()
plt.figure()
(pp, ss) = stats.testCompletedInvertedCumulatives((iHH[0], iSS[0]), method = "AD", plot = True);







ih = ihalop[lblpointsH==672];
matpp.hist(i[:,3], bins=1000, histtype="step", linewidth=3);
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
halopgraph = np.vstack([sortedi, ii]).T
matpp.step(sortedi,ii);
np.savetxt('/home/mtllab/Documents/Haloperidol/1strun/CumulativeI_filtered_all_haloperidol.csv', halopgraph, delimiter=',', newline='\n', fmt='%.5e');


isal = isaline[lblpointsS==672];
matpp.hist(i[:,3], bins=1000, histtype="step");
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
salinegraph = np.vstack([sortedi, ii]).T
matpp.step(sortedi,ii);
np.savetxt('/home/mtllab/Documents/Haloperidol/1strun/CumulativeI_filtered_all_saline.csv', salinegraph, delimiter=',', newline='\n', fmt='%.5e');


ih = ihalop[lblpointsH==1031];
matpp.hist(i[:,3], bins=50, histtype="step");
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
halopgraph = np.vstack([sortedi, ii]).T
matpp.step(sortedi,ii);
np.savetxt('/home/mtllab/Documents/Haloperidol/1strun/CumulativeI_filtered_GPI_haloperidol.csv', halopgraph, delimiter=',', newline='\n', fmt='%.5e');


isal = isaline[lblpointsS==1031];
matpp.hist(i[:,3], bins=50, histtype="step");
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
salinegraph = np.vstack([sortedi, ii]).T
matpp.step(sortedi,ii);
np.savetxt('/home/mtllab/Documents/Haloperidol/1strun/CumulativeI_filtered_GPI_saline.csv', salinegraph, delimiter=',', newline='\n', fmt='%.5e');



ih = ihalop[lblpointsH==536];
#matpp.hist(i[:,3], bins=50, histtype="step");
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
halopgraph = np.vstack([sortedi, ii]).T
matpp.step(sortedi,ii);
#np.savetxt('/home/mtllab/Documents/Haloperidol/1strun/CumulativeI_filtered_GPI_haloperidol.csv', halopgraph, delimiter=',', newline='\n', fmt='%.5e');


isal = isaline[lblpointsS==536];
#matpp.hist(i[:,3], bins=50, histtype="step");
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
salinegraph = np.vstack([sortedi, ii]).T
matpp.step(sortedi,ii);
#np.savetxt('/home/mtllab/Documents/Haloperidol/1strun/CumulativeI_filtered_GPI_saline.csv', salinegraph, delimiter=',', newline='\n', fmt='%.5e');







i = ihalop[lblpointsH==749];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);


i = isaline[lblpointsS==749];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);




i = ihalop[lblpointsH==961];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);


i = isaline[lblpointsS==961];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);



########################################
########################################
import iDISCO.IO.IO as io;
import numpy as np;
import matplotlib.pyplot as matpp;
import iDISCO.Analysis.Label as lbl;
import os, numpy, math

annotationFile = '/home/mtllab/Documents/warping/annotation_25_right.tif';



phalop68, ihalop68 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1268/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1268/intensities.npy')));
phalop69, ihalop69 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1269/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1269/intensities.npy')));
phalop70, ihalop70 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1270/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1270/intensities.npy')));
phalop66, ihalop66 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1266/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1266/intensities.npy')));
phalop67, ihalop67 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1267/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1267/intensities.npy')));


lblpointsH68 = lbl.labelPoints(phalop68, annotationFile);
lblpointsH69 = lbl.labelPoints(phalop69, annotationFile);
lblpointsH70 = lbl.labelPoints(phalop70, annotationFile);
lblpointsH66 = lbl.labelPoints(phalop66, annotationFile);
lblpointsH67 = lbl.labelPoints(phalop67, annotationFile);

phalop71, ihalop71 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1271/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1271/intensities.npy')));
phalop72, ihalop72 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1272/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1272/intensities.npy')));
phalop73, ihalop73 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1273/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1273/intensities.npy')));
phalop74, ihalop74 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1274/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1274/intensities.npy')));
phalop75, ihalop75 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1275/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1275/intensities.npy')));


lblpointsH71 = lbl.labelPoints(phalop71, annotationFile);
lblpointsH72 = lbl.labelPoints(phalop72, annotationFile);
lblpointsH73 = lbl.labelPoints(phalop73, annotationFile);
lblpointsH74 = lbl.labelPoints(phalop74, annotationFile);
lblpointsH75 = lbl.labelPoints(phalop75, annotationFile);





region = 64;
regionName = 'cellsize_AD.csv';


i = ihalop66[lblpointsH66==region];
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii, color='red');
graph = np.vstack([sortedi, ii]).T
np.savetxt(os.path.join('/home/mtllab/Documents/Haloperidol/1266/', regionName), graph, delimiter=',', newline='\n', fmt='%.5e');


i = ihalop67[lblpointsH67==region];
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='red');
graph = np.vstack([sortedi, ii]).T
np.savetxt(os.path.join('/home/mtllab/Documents/Haloperidol/1267/', regionName), graph, delimiter=',', newline='\n', fmt='%.5e');


i = ihalop68[lblpointsH68==region];
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='red');
graph = np.vstack([sortedi, ii]).T
np.savetxt(os.path.join('/home/mtllab/Documents/Haloperidol/1268/', regionName), graph, delimiter=',', newline='\n', fmt='%.5e');

i = ihalop69[lblpointsH69==region];
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='red');
graph = np.vstack([sortedi, ii]).T
np.savetxt(os.path.join('/home/mtllab/Documents/Haloperidol/1269/', regionName), graph, delimiter=',', newline='\n', fmt='%.5e');


#i = ihalop70[lblpointsH70==region];
#sortedi = np.sort(i[:,3]);
#ii = np.arange(sortedi.size);
#matpp.step(sortedi,ii,color='orange');


#i = ihalop71[lblpointsH71==region];
#sortedi = np.sort(i[:,3]);
#ii = np.arange(sortedi.size);
#matpp.step(sortedi,ii,color='blue');
#graph = np.vstack([sortedi, ii]).T
#np.savetxt(os.path.join('/home/mtllab/Documents/Haloperidol/1271/', regionName), graph, delimiter=',', newline='\n', fmt='%.5e');

i = ihalop72[lblpointsH72==region];
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='blue', linewidth = 3);
graph = np.vstack([sortedi, ii]).T
np.savetxt(os.path.join('/home/mtllab/Documents/Haloperidol/1272/', regionName), graph, delimiter=',', newline='\n', fmt='%.5e');

i = ihalop73[lblpointsH73==region];
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='blue');
graph = np.vstack([sortedi, ii]).T
np.savetxt(os.path.join('/home/mtllab/Documents/Haloperidol/1273/', regionName), graph, delimiter=',', newline='\n', fmt='%.5e');

i = ihalop74[lblpointsH74==region];
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='blue');
graph = np.vstack([sortedi, ii]).T
np.savetxt(os.path.join('/home/mtllab/Documents/Haloperidol/1274/', regionName), graph, delimiter=',', newline='\n', fmt='%.5e');

i = ihalop75[lblpointsH75==region];
sortedi = np.sort(i[:,3]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='blue');
graph = np.vstack([sortedi, ii]).T
np.savetxt(os.path.join('/home/mtllab/Documents/Haloperidol/1275/', regionName), graph, delimiter=',', newline='\n', fmt='%.5e');








psaline, isaline = io.readPoints((('/home/mtllab/Documents/Haloperidol/1273/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1273/intensities.npy')));

lblpointsS = lbl.labelPoints(psaline, annotationFile);

i = isaline[lblpointsS==672];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);




i = ihalop[lblpointsH==1031];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);


i = isaline[lblpointsS==1031];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);






i = ihalop[lblpointsH==749];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);


i = isaline[lblpointsS==749];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);




i = ihalop[lblpointsH==961];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);


i = isaline[lblpointsS==961];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);




########################################
########################################
import iDISCO.IO.IO as io;
import numpy as np;
import matplotlib.pyplot as matpp;
import iDISCO.Analysis.Label as lbl;


phalop68, ihalop68 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1268/cells-allpoints.npy'),  ('/home/mtllab/Documents/Haloperidol/1268/intensities-allpoints.npy')));
phalop69, ihalop69 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1269/cells-allpoints.npy'),  ('/home/mtllab/Documents/Haloperidol/1269/intensities-allpoints.npy')));
phalop70, ihalop70 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1270/cells-allpoints.npy'),  ('/home/mtllab/Documents/Haloperidol/1270/intensities-allpoints.npy')));
phalop66, ihalop66 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1266/cells-allpoints.npy'),  ('/home/mtllab/Documents/Haloperidol/1266/intensities-allpoints.npy')));
phalop67, ihalop67 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1267/cells-allpoints.npy'),  ('/home/mtllab/Documents/Haloperidol/1267/intensities-allpoints.npy')));


phalop71, ihalop71 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1271/cells-allpoints.npy'),  ('/home/mtllab/Documents/Haloperidol/1271/intensities-allpoints.npy')));
phalop72, ihalop72 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1272/cells-allpoints.npy'),  ('/home/mtllab/Documents/Haloperidol/1272/intensities-allpoints.npy')));
phalop73, ihalop73 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1273/cells-allpoints.npy'),  ('/home/mtllab/Documents/Haloperidol/1273/intensities-allpoints.npy')));
phalop74, ihalop74 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1274/cells-allpoints.npy'),  ('/home/mtllab/Documents/Haloperidol/1274/intensities-allpoints.npy')));
phalop75, ihalop75 = io.readPoints((('/home/mtllab/Documents/Haloperidol/1275/cells-allpoints.npy'),  ('/home/mtllab/Documents/Haloperidol/1275/intensities-allpoints.npy')));

############################################
############################################



BaseDirectory = '/home/mtllab/Documents/whiskers/exploration/1/oldcounts';
points, intensities = io.readPoints((os.path.join(BaseDirectory, 'cells-allpoints.npy'),  os.path.join(BaseDirectory,  'intensities-allpoints.npy')));
points, intensities = thresholdPoints(points, intensities, threshold = (50, 65000), row = (1,0));
io.writePoints((os.path.join(BaseDirectory, 'cells-50.npy'), os.path.join(BaseDirectory,  'intensities-50.npy')), (points, intensities));


dataSource =  '/home/mtllab/Documents/whiskers/exploration/1/150709_0_8X-cfos_18-29-28/18-29-28_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif'
pointSource= os.path.join(BaseDirectory, 'cells-50.npy');
x = (1200,1700);
y = (1100,1900);
z = (300,800);


data = plt.overlayPoints(dataSource, pointSource, x = x, y = y, z = z, pointColor = None);
io.writeData(os.path.join(BaseDirectory, 'cells-50.tif'), data)

matpp.hist(intensities[:,1], histtype = 'step', bins = 500, range = (10,150))

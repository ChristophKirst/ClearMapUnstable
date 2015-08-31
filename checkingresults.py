# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 16:03:57 2015

@author: mtllab
"""

import iDISCO.Visualization.Plot as plt;
import iDISCO.IO.IO as io;

dataSource =  '/home/mtllab/Documents/Haloperidol/1strun/Haloperidol/150630_0_8X-cfos_18-32-28/18-32-28_0_8X-cfos_UltraII_C00_xyz-Table Z\d{4}.ome.tif'
pointSource= '/home/mtllab/Documents/Haloperidol/1strun/Haloperidol/cells.npy';
x = (1000,1500);
y = (700,1200);
z = (580,1080);


data = plt.overlayPoints(dataSource, pointSource, x = x, y = y, z = z, pointColor = None);
io.writeData('/home/mtllab/Documents/Haloperidol/1strun/Haloperidol/checkpoints_20_filteredI.tif', data)


########################################
########################################
import iDISCO.IO.IO as io;
import numpy as np;
import matplotlib.pyplot as matpp;
import iDISCO.Analysis.Label as lbl;

annotationFile = '/home/mtllab/Documents/warping/annotation_25_right.tif';

phalop, ihalop = io.readPoints((('/home/mtllab/Documents/Haloperidol/1strun/Haloperidol/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1strun/Haloperidol/intensities.npy')));

lblpointsH = lbl.labelPoints(phalop, annotationFile);

psaline, isaline = io.readPoints((('/home/mtllab/Documents/Haloperidol/1strun/saline/cells_transformed_to_Atlas.npy'),  ('/home/mtllab/Documents/Haloperidol/1strun/saline/intensities.npy')));

lblpointsS = lbl.labelPoints(psaline, annotationFile);



i = ihalop[lblpointsH==672];
matpp.hist(i[:,0], bins=1000, histtype="step");


sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii);


i = isaline[lblpointsS==672];
matpp.hist(i[:,0], bins=1000, histtype="step");


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





region = 672;

i = ihalop66[lblpointsH66==region];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii, color='red');

i = ihalop67[lblpointsH67==region];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='red');

i = ihalop68[lblpointsH68==region];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='red');

i = ihalop69[lblpointsH69==region];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='red');

i = ihalop70[lblpointsH70==region];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='orange');


i = ihalop71[lblpointsH71==region];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='blue');

i = ihalop72[lblpointsH72==region];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='blue');

i = ihalop73[lblpointsH73==region];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='blue');

i = ihalop74[lblpointsH74==region];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='blue');

i = ihalop75[lblpointsH75==region];
sortedi = np.sort(i[:,0]);
ii = np.arange(sortedi.size);
matpp.step(sortedi,ii,color='blue');








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


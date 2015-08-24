# -*- coding: utf-8 -*-
"""
Inteface to Illastik pixel classification

Note: 
    - genereate a classifier from illastik 0.5 
    - press 'Train and classify' button and then save as file
    - use file for routines here to classify
    - try to avoid too many features in classifier as classification gets very memory intensive otherwise

Created on Thu Jun  4 14:37:06 2015

@author: ckirst

based on ilastik interface from cell profiler
"""

import sys
import numpy
import h5py
import math

import iDISCO.Settings as settings

def initializeIlastik(path = None):
    """Set system path for illastik installation"""
    if path is None:
        path = settings.IlastikPath;
    sys.path.insert(1, path);

initializeIlastik();

Initialized = False;

#initialize ilastik
try:
    #sys.stdout = sys.stderr = open(os.devnull, "w")
    from ilastik.core.dataMgr import DataMgr, DataItemImage
    from ilastik.modules.classification.core.featureMgr import FeatureMgr
#    from ilastik.modules.classification.core.classificationMgr import ClassificationMgr
    from ilastik.modules.classification.core.features.featureBase import FeatureBase
    from ilastik.modules.classification.core.classifiers.classifierRandomForest import ClassifierRandomForest
    from ilastik.modules.classification.core.classificationMgr import ClassifierPredictThread
    from ilastik.core.volume import DataAccessor
    
    Initialized = True;
    #sys.stdout = old_stdout
except ImportError, ilastikImport:
    #sys.stdout = old_stdout
    print """Ilastik import: failed to import the ilastik libraries, set path in Settings.py!"""
    #raise ilastikImport
    
    
from iDISCO.ImageProcessing.BackgroundRemoval import removeBackground

from iDISCO.ImageProcessing.MaximaDetection import findCenterOfMaxima, findIntensity

from iDISCO.Utils.Timer import Timer
from iDISCO.Utils.ParameterTools import writeParameter

from iDISCO.Visualization.Plot import plotTiling



class IlastikClassifier():
    """Ilastik classifier class to handle ilastik runs using a random forest classifier"""    
    
    def __init__(self):
        
        #self.image_name = '/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/hESCells_DAPI.tif'
        #self.classifier_name = '/home/ckirst/Desktop/' + 'classifier.h5' 
        
        self.classifiers = [];
        self.features = [];        
    
    def loadRandomforest(self, filename = None):
        #if filename is None:
        #   filename = self.classifier_name
        
        prefix = 'classifiers'
        hf = h5py.File(filename,'r')
        if prefix in hf:
            cids = hf[prefix].keys()
            hf.close()
            del hf
            
            classifiers = []
            for cid in cids:
                classifiers.append(ClassifierRandomForest.loadRFfromFile(filename, str(prefix + '/' + cid)))   
            self.classifiers = classifiers
        else:
            raise ImportError('No Classifiers in prefix')
    
    
    def loadFeatures(self, filename = None):
        
        featureItems = []
        hf = h5py.File(filename,'r')
        for fgrp in hf['features'].values():
            featureItems.append(FeatureBase.deserialize(fgrp))
        hf.close()
        del hf
        self.features = featureItems
    
               
    def loadClassifier(self, filename = None):
        """Load a classifier """
        
        self.loadRandomforest(filename)
        self.loadFeatures(filename)
    
    def run(self, image):
        """Run ilastik classifier"""
        
        # Transform input image to ilastik conventions
        # 3d = (time,x,y,z,channel) 
        # 2d = (time,1,x,y,channel)
        ndim = len(image.shape);
        if ndim == 2:
            image.shape = (1,1) + image.shape
        elif len(image.shape) == 3:
            image.shape = (1,) + image.shape
        else:
            raise RuntimeError("Error: set image: failed: image must be 2 or 3d.");
        
        #add singelton dim for channel
        image.shape = image.shape + (1,)
        
        # Create ilastik dataMgr
        dataMgr = DataMgr()
        
        di = DataItemImage('')
        di.setDataVol(DataAccessor(image))
        dataMgr.append(di, alreadyLoaded=True)
        dataMgr.module["Classification"]["classificationMgr"].classifiers = self.classifiers
        
          # Create FeatureMgr
        fm = FeatureMgr(dataMgr, self.features)
        fm.prepareCompute(dataMgr)
        fm.triggerCompute()
        fm.joinCompute(dataMgr)
      
        # Predict with loaded classifier
        classificationPredict = ClassifierPredictThread(dataMgr)
        classificationPredict.start()
        classificationPredict.wait()
        
        if ndim == 2:
            image.shape = image.shape[2:-1];
        else:
            image.shape = image.shape[1:-1];
        
        
        #del dataMgr
        pred = classificationPredict._prediction[0];
        if ndim == 2:
            return pred[0,0,:,:,:]
        else:
            return pred[0,:,:,:,:]


def checkIlastikInitialized():
    global Initialized;
    if not Initialized:
        raise RuntimeError('Ilastik is not initialized, set path in Settings.py');
    



def rescaleToIlastik(img, rescale = None, verbose = False, out = sys.stdout, **parameter):
    # rescale to fit with ilasstik representation of intensities
    timer = Timer();
    writeParameter(out = out, head = 'Rescaling:', rescale = rescale);
    
    if not rescale is None:
        img = img.astype('float32') *  rescale;
        img[img > math.pow(2,8)-1] = math.pow(2,8)-1;
        img = img.astype('uint8');
    else:
        img = img.astype('uint8');
    
    if verbose:
        plotTiling(img);
    
    out.write(timer.elapsedTime(head = 'Rescaling') + '\n');  
    return img


def detectCells(img, classifier = None, verbose = False, out = sys.stdout, **parameter):
    """Detect Cells Using a trained classifier in Ilastik"""

    checkIlastikInitialized();
    
    # normalize data
    img2 = rescaleToIlastik(img, verbose = verbose, out = out, **parameter);
    
    #remove backgroudn
    img2 = removeBackground(img, verbose = verbose, out = out, **parameter);
      

    #classify image / assume class 1 are the cells !  
    timer = Timer();  
    writeParameter(out = out, head = 'Ilastik classification:', classifier = classifier);
    
    cls = IlastikClassifier();
    cls.loadClassifier(classifier);
    imgmax = cls.run(img2);
    #print imgmax.shape
    #max probability gives final class
    imgmax = numpy.argmax(imgmax, axis = -1);
    
    # class 1 is used as cells 
    imgmax = imgmax == 1; # class 1 is used as cells 
    
    if verbose:
        plotTiling(imgmax * 1000);
    
    out.write(timer.elapsedTime(head = 'Ilastik classification') + '\n');    
    
    
    #center of maxima
    centers = findCenterOfMaxima(img2, imgmax, verbose = verbose, out = out);
    
    
    #intensity of cells
    cintensity = findIntensity(img, centers, verbose = verbose, out = out, **parameter);

    #intensity of cells in filtered image
    cintensity2 = findIntensity(img2, centers, verbose = verbose, out = out, **parameter);

    
    return ( centers, numpy.vstack((cintensity, cintensity2)).transpose());    
    




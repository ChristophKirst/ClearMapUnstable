# -*- coding: utf-8 -*-
"""
Inteface to Illastik pixel classification

Note: genereate a classifier from illastik 0.5 
      Press 'Train and classify' button and save as file, use file name in routine to classify

      Ilastik classification is very memory intensive !

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
import numpy
import h5py
import math


try:
    #sys.stdout = sys.stderr = open(os.devnull, "w")
    from ilastik.core.dataMgr import DataMgr, DataItemImage
    from ilastik.modules.classification.core.featureMgr import FeatureMgr
#    from ilastik.modules.classification.core.classificationMgr import ClassificationMgr
    from ilastik.modules.classification.core.features.featureBase import FeatureBase
    from ilastik.modules.classification.core.classifiers.classifierRandomForest import ClassifierRandomForest
    from ilastik.modules.classification.core.classificationMgr import ClassifierPredictThread
    from ilastik.core.volume import DataAccessor
    #sys.stdout = old_stdout
except ImportError, ilastikImport:
    #sys.stdout = old_stdout
    print """ilastik import: failed to import the ilastik library. Please follow the instructions on "http://www.ilastik.org" to install ilastik"""
    raise ilastikImport

import cv2

from scipy.ndimage.measurements import label, center_of_mass

from iDISCO.ImageProcessing.Filter.StructureElement import structureElement
from iDISCO.Parameter import ImageProcessingParameter

from iDISCO.Utils.Timer import Timer
from iDISCO.Visualization.Plot import plotTiling, plotOverlayLabel

class IlastikClassifier():
    
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


def rescaleToIlastik(img, parameter = ImageProcessingParameter(), verbose = False):
    # rescale to fit with ilasstik representation
    if parameter.Parameter.Rescale != None:
        img = img.astype('float32') *  parameter.Parameter.Rescale;
        img[img > math.pow(2,8)-1] = math.pow(2,8)-1;
        img = img.astype('uint8');
    else:
        img = img.astype('uint8');
    
    if verbose:
        plotTiling(img);
    
    return img


def detectCells(img, verbose = False, out = sys.stdout, parameter = ImageProcessingParameter()):
    """Detect Cells Using a trained classifier in Ilastik"""

    timer = Timer();
    
    # normalize data
    img = rescaleToIlastik(img, parameter = parameter, verbose = verbose);
    
    if parameter.Parameter.Background != None:     # background subtraction in each slice
        timer.reset();
        se = structureElement('Disk', parameter.Parameter.Background).astype('uint8');
        for z in range(img.shape[2]):
             img[:,:,z] = img[:,:,z] - cv2.morphologyEx(img[:,:,z], cv2.MORPH_OPEN, se)    
        out.write(timer.elapsedTime(head = 'Background') + '\n');
    
        if verbose:
            plotTiling(10*img);       
    
    
    #classify image / assume class 1 are the cells !    
    timer.reset();
    cls = IlastikClassifier();
    cls.loadClassifier(parameter.Parameter.Classifier);
    imgmax = cls.run(img);
    #print imgmax.shape
    #max probability gives final class
    imgmax = numpy.argmax(imgmax, axis = -1);
    
    # class 1 is used as cells 
    imgmax = imgmax == 1; # class 1 is used as cells 
    out.write(timer.elapsedTime(head = 'Ilastik classification') + '\n');
    
    if verbose:
        plotTiling(imgmax * 1000);
    
    #center of maxima
    timer.reset();
    imglab, nlab = label(imgmax);
    
    if nlab > 0:
        centers = numpy.array(center_of_mass(img, imglab, index = numpy.arange(1, nlab)));    
        out.write(timer.elapsedTime(head = 'Cell Centers') + '\n');
    
        if verbose:        
            #plotOverlayLabel(img * 0.01, imglab, alpha = False);
            #plotTiling(img)
            imgc = numpy.zeros(img.shape);
            for i in range(centers.shape[0]):
                imgc[centers[i,0], centers[i,1], centers[i,2]] = 1;
            plotOverlayLabel(img * 0.01, imgc, alpha = False);
            #plotOverlayLabel(img, imgmax.astype('int64'), alpha = True)     
    
        #return centers, imglab, mask
        cintensity = numpy.array([img[centers[i,0], centers[i,1], centers[i,2]] for i in range(centers.shape[0])]);        
    
        return ( centers, cintensity );
        
    else:
        out.write(timer.elapsedTime(head = 'Cell Centers') + '\n');
        out.write('No Cells found !');
        return ( numpy.zeros((0,3)), numpy.zerso(0) );    
    
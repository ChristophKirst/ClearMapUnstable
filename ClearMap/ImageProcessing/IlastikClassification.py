# -*- coding: utf-8 -*-
"""
Inteface to Illastik pixel classification

This module allows to integrate ilastik pixel classification into the *ClearMap*
pipeline. 

To train a classifier ilastik 0.5 should be used following these steps:

    * genereate a classifier from illastik 0.5 
    * press 'Train and classify' button (eventhough the online traning is running) !
    * save the classifier to a file
    * use the classifiers file name in the ClearMap routine :func:`classifyPixel`
    * try to avoid too many features in the classifier 
      as classification gets very memory intensive otherwise
    
Note:
    Note that ilastik classification works in parallel, thus it is advised to 
    process the data sequentially, see 
    :fun:`~ClearMap.Imageprocessing.StackProcessing.sequentiallyProcessStack`  

Note:
    Ilastik 0.5 works for images in uint8 format !

References:
    `Ilastik <http://ilastik.org/>`_
    
Author
""""""
   - Based on the ilastik interface from `cell profiler <http://www.cellprofiler.org/>`_
   - Christoph Kirst, The Rockefeller University, New York City, 2015
"""

import sys
import numpy
import h5py
import math

import ClearMap.Settings as settings

#from ClearMap.ImageProcessing.BackgroundRemoval import removeBackground
from ClearMap.ImageProcessing.MaximaDetection import findCenterOfMaxima, findIntensity
from ClearMap.ImageProcessing.StackProcessing import writeSubStack

from ClearMap.Utils.Timer import Timer
from ClearMap.Utils.ParameterTools import getParameter, writeParameter

from ClearMap.Visualization.Plot import plotTiling


def initializeIlastik(path = None):
    """Set system path for illastik installation"""
    if path is None:
        path = settings.IlastikPath;
    sys.path.insert(1, path);

initializeIlastik();

Initialized = False;
"""bool: True if ilastik interface was sucessfully initialized ."""


#initialize ilastik
try:
    #sys.stdout = sys.stderr = open(os.devnull, "w")
    from ilastik.core.dataMgr import DataMgr, DataItemImage
    from ilastik.modules.classification.core.featureMgr import FeatureMgr
    #from ilastik.modules.classification.core.classificationMgr import ClassificationMgr
    from ilastik.modules.classification.core.features.featureBase import FeatureBase
    from ilastik.modules.classification.core.classifiers.classifierRandomForest import ClassifierRandomForest
    from ilastik.modules.classification.core.classificationMgr import ClassifierPredictThread
    from ilastik.core.volume import DataAccessor
    
    Initialized = True;
    #sys.stdout = old_stdout
    
except ImportError, ilastikImport:
    #sys.stdout = old_stdout
    print """Ilastik import: failed to import the ilastik libraries, set path in ClraMap.Settings!"""
    #raise ilastikImport
    
    

class _IlastikClassifier():
    """Class to handle ilastik pixel classification

    Attributes:
        classifiers 
        features
    """    
    
    def __init__(self):
        
        #self.image_name = '/home/ckirst/Science/Simulation/Matlab/StemCell3D/Test/Images/hESCells_DAPI.tif'
        #self.classifier_name = '/home/ckirst/Desktop/' + 'classifier.h5' 
        
        self.classifiers = [];
        self.features = [];        
    
    def loadRandomforest(self, filename = None):
        """Load a random forest classifier
        
        Arguments:
            filename (str): ilastik classifier file
        """
        
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
        """Load features
        
        Arguments:
            filename (str): ilastik classifier file
        """
        
        featureItems = []
        hf = h5py.File(filename,'r')
        for fgrp in hf['features'].values():
            featureItems.append(FeatureBase.deserialize(fgrp))
        hf.close()
        del hf
        self.features = featureItems
    
               
    def loadClassifier(self, filename = None):
        """Load a classifier 
        
        Arguments:
            filename (str): ilastik classifier file
        """
        
        self.loadRandomforest(filename)
        self.loadFeatures(filename)
        
    
    def run(self, image):
        """Run ilastik classifier
        
        Arguments:
            image (array):  image data
            
        Returns:
            array: probabilities of belonging to a class 
            (image.shape, number of classes)
        """
        
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


def _checkIlastikInitialized():
    """Check if Ilastik is initialized"""
    global Initialized;
    if not Initialized:
        raise RuntimeError('Ilastik is not initialized, set path in Settings.py');
  


  
def rescaleToIlastik(img, rescale = None, verbose = False, out = sys.stdout, **parameter):
    """Rescale image to achieve uint8 format used by ilasstik 
    
    The function rescales the image and converts the image to uin8
    to fit with the image representation used by ilastik.
    
    Arguments:
        img (array): image data
        rescale (float or None): rescaling factor
    
    Returns:
        array : uint8 version of the image 
    """
    
    if not rescale is None:
        img = img.astype('float32') *  rescale;
        img[img > math.pow(2,8)-1] = math.pow(2,8)-1;
        img = img.astype('uint8');
    else:
        img = img.astype('uint8');
    
    return img



def classifyPixel(img, classifyPixelParameter = None, subStack = None, verbose = False, out = sys.stdout, **parameter):
    """Detect Cells Using a trained classifier in Ilastik
    
    Arguments:
        img (array): image data
        classifyPixelParameter (dict):
            ============ ==================== ===========================================================
            Name         Type                 Descritption
            ============ ==================== ===========================================================
            *classifier* (str or  None)       saves result of labeling the differnet maxima
                                              if None dont correct image for illumination, if True the 
            *rescale*    (float or None)      optional rescaling of the image to fit uint8 format 
                                              used by ilastik
            *save*       (str or None)        save the propabilities to belong to the classes to a file
            *verbose*    (bool or int)        print / plot information about this step 
            ============ ==================== ===========================================================
        subStack (dict or None): sub-stack information 
        verbose (bool): print progress info 
        out (object): object to write progress info to
    
    Returns:
        array: probabilities for each pixel to belong to a class in the 
               classifier, shape is (img.shape, number of classes)
    """

    _checkIlastikInitialized();
    
    classifier = getParameter(classifyPixelParameter, "classifier", None);
    rescale    = getParameter(classifyPixelParameter, "rescale", None);    
    save       = getParameter(classifyPixelParameter, "save", None);   
     
    if verbose:
        writeParameter(out = out, head = 'Ilastik classification:', classifier = classifier, rescale = rescale, save = save);        
    
    
    timer = Timer(); 
    
    # normalize data
    img2 = rescaleToIlastik(img, verbose = verbose, out = out, **parameter);
    
    #remove background
    #img2 = removeBackground(img, verbose = verbose, out = out, **parameter);
      
    #classify image

    if classifier is None:        
        return img2;
    
    cls = _IlastikClassifier();
    cls.loadClassifier(classifier);
    imgclass = cls.run(img2);
    
    if not save is None:
        for i in range(imgclass.shape[4]):
            fn = save[:-4] + '_class_' + str(i) + save[-4:];
            writeSubStack(fn, imgclass[:,:,:,i], subStack = subStack)
      
    if verbose > 1:
        for i in range(imgclass.shape[4]):
            plotTiling(imgclass[:,:,:,i]);
    
    if verbose:
        out.write(timer.elapsedTime(head = 'Ilastik classification') + '\n');    
    
    return imgclass;


def classifyCells(img, classifyCellsParameter = None, classifier = None, rescale = None, save = None, verbose = False,
                  subStack = None, out = sys.stdout, **parameter):
    """Detect Cells Using a trained classifier in Ilastik
    
    The routine assumes that the first class is identifying the cells.
        
    Arguments:    
       img (array): image data
        classifyPixelParameter (dict):
            ============ ==================== ===========================================================
            Name         Type                 Descritption
            ============ ==================== ===========================================================
            *classifier* (str or  None)       saves result of labeling the differnet maxima
                                              if None dont correct image for illumination, if True the 
            *rescale*    (float or None)      optional rescaling of the image to fit uint8 format 
                                              used by ilastik
            *save*       (str or None)        save the detected cell pixel to a file
            *verbose*    (bool or int)        print / plot information about this step 
            ============ ==================== ===========================================================
        subStack (dict or None): sub-stack information 
        verbose (bool): print progress info 
        out (object): object to write progress info to
    
    Returns:
        tuple: centers of the cells, intensity measurments
        
    Note:    
        The routine could be poteNtially refined to make use of background 
        detected by ilastik
    """

    _checkIlastikInitialized();
    
    classifier = getParameter(classifyCellsParameter, "classifier", classifier);
    rescale    = getParameter(classifyCellsParameter, "rescale", rescale);    
    save       = getParameter(classifyCellsParameter, "save", save);   
    verbose    = getParameter(classifyCellsParameter, "verbose", verbose);
     
    if verbose:
        writeParameter(out = out, head = 'Ilastik cell detection:', classifier = classifier, rescale = rescale, save = save);        

    timer = Timer(); 

    _checkIlastikInitialized();
    
    # normalize data
    img2 = rescaleToIlastik(img, verbose = verbose, out = out, **parameter);
    
    #remove background
    #img2 = removeBackground(img, verbose = verbose, out = out, **parameter);
      
    #classify image / assume class 1 are the cells !  
    timer = Timer();  
    
    cls = _IlastikClassifier();
    cls.loadClassifier(classifier);
    imgmax = cls.run(img2);
    #print imgmax.shape
    #max probability gives final class
    imgmax = numpy.argmax(imgmax, axis = -1);
    
    # class 1 is used as cells 
    imgmax = imgmax == 1; # class 1 is used as cells 
    
    if save:
        writeSubStack(save, imgmax, subStack = subStack)
      
    if verbose > 1:
        plotTiling(imgmax);
        
    #center of maxima
    centers = findCenterOfMaxima(img2, imgmax, verbose = verbose, out = out, **parameter);
    
    #intensity of cells
    cintensity = findIntensity(img, centers, verbose = verbose, out = out, **parameter);

    #intensity of cells in filtered image
    cintensity2 = findIntensity(img2, centers, verbose = verbose, out = out, **parameter);
    
    if verbose:
        out.write(timer.elapsedTime(head = 'Ilastik cell detection') + '\n');    
    
    return ( centers, numpy.vstack((cintensity, cintensity2)).transpose() );    
    



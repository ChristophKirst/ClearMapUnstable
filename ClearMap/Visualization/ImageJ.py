# -*- coding: utf-8 -*-
"""
View images in ImageJ

Supported functionality:

    * open image stack in imagej
    * opem image stack in 3DViewer
    
Note: 
  This is a simple interface operating via command line arguments to imagej.

"""
#:copyright: Copyright 2015 by Christoph Kirst, The Rockefeller University, New York City
#:license: GNU, see LICENSE.txt for details.

import numpy
import os
import tempfile
import glob

import ClearMap.Settings as settings
import ClearMap.IO as io


##############################################################################
### Initialization and Enviroment Settings
##############################################################################

ImageJBinary = None;
"""str: the ImageJ executable

Notes:
    - setup in :func:`initializeImageJ`
"""
    
Initialized = False;
"""bool: True if the ImageJ binary is setup

Notes:
    - setup in :func:`initializeImageJ`
"""

    
def printSettings():
    """Prints the current ImageJ configuration
    
    See also:
        :const:`ImageJBinary`, :const:`Initialized`
    """
    
    global ImageJBinary, Initialized
    
    if Initialized:
        print "ImageJBinary = %s" % ImageJBinary;
    else:
        print "ImageJ not initialized";



def initializeImageJ(path = None):
    """Initialize settings for the ImageJ

    Arguments:
        path (str or None): path to ImageJ root directory, if None 
        :const:`ClearMap.Settings.ImageJPath` is used.
        
    See also:
        :const:`ImageJBinary`,
        :const:`Initialized`
    """
    
    global ImageJBinary, Initialized
    
    if path is None:
        path = settings.ImageJPath;
    
    
    #search for elastix binary
    imagejbin = os.path.join(path, 'ImageJ-linux64');
    if os.path.exists(imagejbin):
        ImageJBinary = imagejbin;
    else:
        raise RuntimeError("Cannot find ImageJ binary %s, set path in Settings.py accordingly!" % imagejbin);
    
    Initialized = True;
    
    print "ImageJ sucessfully initialized from path: %s" % path;
    
    return path;



initializeImageJ();


def checkImageJInitialized():
    """Checks if ImageJ is initialized
    
    Returns:
        bool: True if ImageJ paths are set.
    """
    
    global Initialized;
    
    if not Initialized:
        raise RuntimeError("ImageJ not initialized: run initializeImageJ(path) with proper path to ImageJ first");
    #print ElastixSettings.ElastixBinary;

    return True;


##############################################################################
### Run ImageJ
##############################################################################

def openData(dataSource, x = all, y = all, z = all, inverse = False, cleanUp = True): 
    """Open image in ImageJ 
    
    Arguments:
        dataSouce (str or array): volumetric image data
        x, y, z (all or tuple): sub-range specification
        inverse (bool):invert image
    
    Returns:
        (object): figure handle
    """

    checkImageJInitialized();
    
    if isinstance(dataSource, numpy.ndarray):
      filename = tempfile.mktemp(suffix = '.mhd', prefix = 'CM_ImageJ');
      io.writeData(filename, dataSource, x = x, y = y, z = z);
      temp= True
      dSize = dataSource.shape;
    else:
      filename = dataSource;
      temp = False;
      dSize = io.dataSize(dataSource);
    
    colorImage = len(dSize) == 4;
  
    if colorImage:
      macro = ('open("%s"); ' % filename) + \
               'run("Stack to Hyperstack...", "order=xyzct channels=%d slices=%d frames=1 display=Color"); ' % (dSize[3], dSize[2]) + \
               'Stack.setDisplayMode("composite");';
    else:
      macro = ('open("%s");' % filename);
    
    cmd = ImageJBinary + " -eval '%s'" % macro;
    
    print 'running: %s' % cmd
    res = os.system(cmd);
  
    if res != 0:
      raise RuntimeError('openData: failed executing: ' + cmd);
    
    
    if cleanUp and temp:
      os.remove(filename);
      filepath, imagename = os.path.split(filename);
      imagename = imagename[:-4] + '.raw';
      os.remove(os.path.join(filepath, imagename));
    
    return macro;


def openData3D(dataSource, x = all, y = all, z = all, cleanUp = True):
    """Open image in ImageJ 
    
    Arguments:
        dataSouce (str or array): volumetric image data
        x, y, z (all or tuple): sub-range specification
        inverse (bool):invert image
    
    Returns:
        (object): figure handle
    """

    checkImageJInitialized();
    
    if isinstance(dataSource, numpy.ndarray):
      filename = tempfile.mktemp(suffix = '.mhd', prefix = 'CM_ImageJ');
      io.writeData(filename, dataSource, x = x, y = y, z = z); 
      filepath, imagename = os.path.split(filename);
      imagename = imagename[:-4] + '.raw';
      temp = True; 
      dSize = dataSource.shape;
    else:
      filename = dataSource;
      filepath, imagename = os.path.split(filename);
      temp = False;
      
    if len(dSize) == 4:
      colorImage = True;
    else:
      colorImage = False;
      
    if colorImage:
      macro = ('open("%s"); ' % filename) + \
               'run("Stack to Hyperstack...", "order=xyzct channels=%d slices=%d frames=1 display=Color"); ' % (dSize[3], dSize[2]) + \
               'Stack.setDisplayMode("composite"); ' + \
               'run("3D Viewer"); call("ij3d.ImageJ3DViewer.setCoordinateSystem", "false"); ' + \
               'call("ij3d.ImageJ3DViewer.add", "%s", "None", "%s", "0", "true", "true", "true", "1", "0");' % (imagename, imagename);
    else:
      macro = ('open("%s");' % filename) + ' run("RGB Color"); run("3D Viewer"); call("ij3d.ImageJ3DViewer.setCoordinateSystem", "false"); ' + \
               'call("ij3d.ImageJ3DViewer.add", "%s", "None", "%s", "0", "true", "true", "true", "1", "0");' % (imagename, imagename);
    
    cmd = ImageJBinary + " -eval '%s'" % macro;
    
    print 'running: %s' % cmd
    res = os.system(cmd);
  
    if res != 0:
      raise RuntimeError('openData3D: failed executing: ' + cmd);
    
    if cleanUp and temp:
      os.remove(filename);
      os.remove(os.path.join(filepath, imagename));
    
    return macro;

    
def cleanTemporary():
  """Cleans temporary folder from CM_ImageJ* files"""
  
  fns = glob.glob(os.path.join(tempfile.tempdir, 'CM_ImageJ*'));
  for f in fns:
    os.remove(f);
 

   
def test():
  
  import ClearMap.Visualization.ImageJ as ij
  import numpy as np
  import glob 
  
  reload(ij)

  data = np.random.rand(20, 20, 20);
  ij.openData(data)

  reload(ij)
  ij.openData3D(data)
  
  print glob.glob(os.path.join(tempfile.tempdir, 'CM_ImageJ*'))
  ij.cleanTemporary();
  print glob.glob(os.path.join(tempfile.tempdir, 'CM_ImageJ*'))


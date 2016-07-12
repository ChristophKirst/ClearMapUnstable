# -*- coding: utf-8 -*-
"""
Interface to TeraStitcher for stitching tiled large volumetric images

The TeraStitcher documentation can be found `here <http://abria.github.io/TeraStitcher/>`_.

Main routine is :func:`stitchData`
    
See Also:
    :mod:`~ClearMap.Alignment`
"""
#:copyright: Copyright 2015 by Christoph Kirst, The Rockefeller University, New York City
#:license: GNU, see LICENSE.txt for details.


import os
import tempfile
import shutil
import numpy as np
import re
import natsort

import ClearMap.Settings as settings
import ClearMap.IO as io

##############################################################################
### Initialization and Enviroment Settings
##############################################################################

TeraStitcherBinary = None;
"""str: the TeraStitcher executable

Notes:
    - setup in :func:`initializeStitcher`
"""
    
Initialized = False;
"""bool: True if the TeraStitcher binary is setup

Notes:
    - setup in :func:`initializeStitcher`
"""

    
def printSettings():
    """Prints the current TeraStitcher configuration
    
    See also:
        :const:`TeraStitcherBinary`, :const:`Initialized`
    """
    
    global TeraStitcherBinary, Initialized
    
    if Initialized:
        print "TeraStitcherBinary = %s" % TeraStitcherBinary;
    else:
        print "TeraStitcher not initialized";



def initializeStitcher(path = None):
    """Initialize settings for the TeraStitcher

    Arguments:
        path (str or None): path to TeraStitcher root directory, if None 
        :const:`ClearMap.Settings.TeraStitcherPath` is used.
        
    See also:
        :const:`TeraStitcherBinary`,
        :const:`Initialized`
    """
    
    global TeraStitcherBinary, Initialized
    
    if path is None:
        path = settings.TeraStitcherPath;
        
    print path
    
    #search for elastix binary
    terasticherbin = os.path.join(path, 'bin/terastitcher');
    if os.path.exists(terasticherbin):
        TeraStitcherBinary = terasticherbin;
    else:
        raise RuntimeError("Cannot find TeraSticher binary %s, set path in Settings.py accordingly!" % terasticherbin);
    
    Initialized = True;
    
    print "TeraSticher sucessfully initialized from path: %s" % path;
    
    return path;



initializeStitcher();


def checkSticherInitialized():
    """Checks if TeraSticher is initialized
    
    Returns:
        bool: True if TeraSticher paths are set.
    """
    
    global Initialized;
    
    if not Initialized:
        raise RuntimeError("TeraSticher not initialized: run initializeTeraSticher(path) with proper path to TeraSticher first");
    #print ElastixSettings.ElastixBinary;

    return True;



##############################################################################
### Basic interface routines
##############################################################################


def findFileList(filename, sort = True, groups = None, absolute = True):
  """Returns the list of files the match the regular expression (including variable paths)
  
  Arguments:
      filename (str): file name as regular expression 
      sort (bool): naturally sort the list
      groups (list or None): if list also return a values found for the indicated group identifier
  
  Returns:
      list: file names
      list: if groups is not None, a list of lists of group identifier
  """
  
  if sort:
    def listdir(path):
      for f in natsort.natsorted(os.listdir(path)):
        if not f.startswith('.'):
          yield f;
  else:
    def listdir(path):
      for f in os.listdir(path):
        if not f.startswith('.'):
          yield f;
  
  # split full paths in case path contains regular expression
  fnsp = filename.split(os.path.sep);

  # handle absolute path 
  if fnsp[0] == '': # aboslute path
    files = [os.path.sep];
    fnsp.pop(0);
  else: # relative path
    if absolute:
      files = [os.path.abspath(os.curdir)];
    else:
      files = ['.'];
  
  # handle pure directory expression -> all files
  if len(fnsp) > 0 and fnsp[-1] == '':
      fnsp[-1] = '.*';
  
  if groups is None:  
    #no group info 
    for n in fnsp:    
      search = re.compile(n).search;
      newfiles = [];
      for f in files:
        matches = map(search, listdir(f));
        matches = [x for x in matches if x];
        for m in matches:
          if f.endswith(os.path.sep):
             newfiles.append(f + m.string);
          else:
            newfiles.append(f + os.path.sep + m.string);
      files = newfiles;
    
    return files;
  
  else: #with group info
    infos = [[None for x in groups]];
    for n in fnsp:    
      search = re.compile(n).search;
      newfiles = [];
      newinfos = [];
      
      for f,i in zip(files, infos):
        matches = map(search, listdir(f));
        matches = [x for x in matches if x];
        for m in matches:
          if f.endswith(os.path.sep):
            newfiles.append(f + m.string);
          else:
            newfiles.append(f + os.path.sep + m.string);
          
          d = m.groupdict();
          ii = list(i);
          for k,v in d.iteritems():
            try:
              j = groups.index(k);
              ii[j] = v;
            except:
              pass
          newinfos.append(ii);
      
      files = newfiles;
      infos = newinfos;
    
    return files, infos;



def findFirstFile(filename, hidden = False, sort = True):
  """Returns the first file that matches the regular expression
  
  Arguments:
      filename (str): file name as regular expression
      hidden (bool): if True include hidden files
      sort (bool): naturally sort the files if True, if False and arbirary file that matched first is returned
      
  Returns:
      string or None: first file name that matches the regular expression or None if no match found
    
  Note:
  
    For large file numbers speed can be imporved by setting the option sort = False
  """
  
  if hidden:
    if sort:
      def listdir(path):
        return natsort.natsorted(os.listdir(path));
    else:
      listdir = os.listdir
  else:
    if sort:
      def listdir(path):
        for f in natsort.natsorted(os.listdir(path)):
          if not f.startswith('.'):
            yield f;
    else:
      def listdir(path):
        for f in os.listdir(path):
          if not f.startswith('.'):
            yield f;
  
  # split full paths in case path contains regular expression
  fnsp = filename.split(os.path.sep);

  # handle absolute path 
  if fnsp[0] == '': # aboslute path
    fn = os.path.sep;
    fnsp.pop(0);
  else: # relative path
    if absolute:
      fn = os.path.abspath(os.curdir);
    else:
      fn = '.';
  
  #no group info 
  for n in fnsp:    
    search = re.compile(n).search;
    match = None;
    for x in listdir(fn):
        if search(x):
          match = x;
          break;  
    
    if match is not None:
        if fn.endswith(os.path.sep):
          fn = fn + match;
        else:
          fn = fn + os.path.sep + match  
    else:
      return None;
  
  return fn;


def findFileInfo(filename):
  """Tries to infer relevant information from filename for tiling / stitching
  
  Arguments:
    filename (str): filename to infer information from
    
  Returns:
    dict: dictionary with relavant information: resolution (microns per pixel), overlap (microns), size (pixel)
    
  Note:
    resolution is in microns per pixel, overlap is in microns and size is in pixel
  """
  #TODO: move this to the individual file readers  
  
  #this is for ome tif images
  from PIL import Image
  from PIL.ExifTags import TAGS
  from lxml import etree as et
   
  def ome_info(fn):
      ret = {}
      i = Image.open(fn)
      info = i.tag.as_dict();
      for tag, value in info.iteritems():
          decoded = TAGS.get(tag)
          ret[decoded] = value
      return ret
      
  imginfo = ome_info(filename);
  keys = imginfo.keys();
  
  finfo = {'resolution' : None, 'overlap' : None, 'size' : None};
  
  # get image sizes  
  if ('ImageHeight' in keys) and ('ImageWidth' in keys):
    finfo['size'] = (imginfo['ImageWidth'], imginfo['ImageHeight']);
  else:
    finfo['size'] = io.dataSize(filename)
  
  
  if 'ImageDescription' in keys:
    imgxml = imginfo['ImageDescription'];
    if isinstance(imgxml, tuple):
      imgxml = imgxml[0];
    imgxml = et.fromstring(str(imgxml));
    
    # get resolution
    pix = [x for x in imgxml.iter('{*}Pixels')];
    if len(pix)> 0:
      pix = pix[0].attrib;
      keys = pix.keys();
      if 'PhysicalSizeX' in keys and 'PhysicalSizeY' in keys and 'PhysicalSizeZ' in keys:
        finfo['resolution'] = (float(pix['PhysicalSizeX']), float(pix['PhysicalSizeY']), float(pix['PhysicalSizeZ']));
  
    # get overlap
    e1 = [x for x in imgxml.iter('{*}xyz-Table_X_Overlap')];
    e2 = [x for x in imgxml.iter('{*}xyz-Table_Y_Overlap')];
    if len(e1) > 0 and len(e2) > 0:
      finfo['overlap'] = (float(e1[0].attrib['Value']), float(e2[0].attrib['Value']));
      
  return finfo;


def displacements(size, resolution = (1.0, 1.0),  overlap = (0.0, 0.0), units = 'Microns'):
  """Calculates the displacements of the tiles in microns given resolution, overlaps and image size
  
  Arguments:
    size (tuple): image size in pixel 
    resolution (tuple): image resolution in microns per pixel
    overlap (tuple): overlap of the images in microns
    units (str): 'Pixel' or 'Microns'
    
  Returns:
    tuple: displacements in x and y of the tiles
  """

  size = np.array(size, dtype = float); 
  size = size[:2];
  resolution = np.array(resolution, dtype = float);
  resolution = resolution[:2];
  overlap = np.array(overlap, dtype = float);
  overlap = overlap[:2];
  
  overlap_pix = overlap / resolution;
  
  if units == 'Pixel':
    return (size - overlap_pix);
  else:    
    return (size - overlap_pix) * resolution;
  
  
def xmlImportStack(row, col, zrange, displacement, stitchable = False, directory = None, expression = None): 
  """Creates a single stack xml information for the xmlf import file of TeraSticher
  
  Arguments:
    row, col (int or str): the row and column specifications
    zrange (tuple or str): the z range specifications
    displacement (tuple): the (x,y,z) displacements of the stack in pixel / voxel
    stitchable (bool or str):  stitachable parameter
    directory (str): directory of stack
    expression (str or None): regular expression specifing the images of the (row,col) tile
  
  Returns:
    etree: xml entry as lxml etree
  """
    
  s = etree.Element('Stack');
  for e in ['NORTH_displacements', 'EAST_displacements', 'SOUTH_displacements', 'WEST_displacements']:
    s.append(etree.Element(e))
  
  s.set('ROW', str(row))
  s.set('COL', str(col))
  s.set('ABS_V', str(displacement[0]));
  s.set('ABS_H', str(displacement[1]));
  s.set('ABS_D', str(displacmeent[2])); 
  s.set('STITCHABLE', str(stitchable).lower());
  if directory is None or directory == '':
    directory = '.';
  s.set('DIR_NAME', str(directory));
  if expression is None:
    expression = '';
  s.set('IMG_REGEX', str(expression));
  if isinstance(zrange, tuple) or isinstance(zrange, list):
    zrange = '[%d, %d)' % (zrange[0], zrange[1]);
  s.set('Z_RANGES', zrange)
  
  return s;


def xmlImport(baseDirectory, size, resolution = (1.0, 1.0, 1.0), 
              origin = (0.0, 0.0, 0.0), overlap = (0.0, 0.0), tiling = (1,1), 
              tileExpression = None, zRange = None, 
              xmlFile = None, asString = False):
  """Creates an xml import file specifing the data for the TeraStitcher
  
  Arguments:
    baseDirectory (str): base directory for the volumetric data
    size (tuple): size of the volume in pixel
    resolution (tuple): resolution of the image n microns per pixel
    origin (tuple): origin of the image
    overlap (tuple): overlap of the tiles in x,y in microns
    tiling (tuple): tile dimensions
    tileExpression (function or None): tileExpression(row,col) should return a 
                                       regular expression specifing the path and images 
                                       of the (row,col) tile, if None all 
                                       images in the tileirectory are used
    zRange (function or None): zRange(row, col) should return the zrange for 
                               the (row,col) tile as a string, if None, full 
                               range is used
    xmlFile (str or None): file name to save the xml to, if None return xml 
                           string or xml element tree
    asString (bool): if True return as string otherwise as xml element tree
  """
  
  xml = etree.Element('TeraStitcher');
  
  e = etree.Element('stacks_dir');
  e.set('value', baseDirectory);  
  xml.append(e);
  
  e = etree.Element('voxel_dims');
  e.set('V', str(resolution[0]));
  e.set('H', str(resolution[1]));
  e.set('D', str(resolution[2]));
  xml.append(e);
  
  e = etree.Element('origin');
  e.set('V', str(origin[0]));
  e.set('H', str(origin[1]));
  e.set('D', str(origin[2]));
  xml.append(e);

  e = etree.Element('mechanical_displacements');
  displ = displacements(size, resolution, overlap, form = 'Microns');
  e.set('V', str(displ[0]));
  e.set('H', str(displ[1]));
  xml.append(e);
  
  e = etree.Element('dimensions');
  e.set('stack_columns', str(tiling[0]));
  e.set('stack_rows', str(tiling[1]));
  e.set('stack_slices', str(size[2]));
  xml.append(e);  
  
  if zRange is all or zRange is None:
    zRange = lambda row,col: (0, size[2]); 
  
  displ = displacements(size, resolution, overlap, form = 'Pixel');
  
  if tileDirectory is None:
    tileDirectory = lambda row,col: '.';
  
  if tileExpression is None:
    tileExpression = lambda row,col: '';
  
  stack = e.etree.Element('STACKS');
  for row in range(tiling[0]):
    for col in range(tiling[1]):
      stackdispl = [displ[0] * row, displ[1] * col];
      te = tileExpression(row,col);
      dn,fn = os.path.split(te);      
      e = xmlImportStack(row, col, zRange(row,col), stackdispl, stitchable = False, directory = dn, expression = fn);
      stack.append(e);
  
  xml.append(stack);
  
  
  if xmlFile is None:
    if asString:
      return etree.tostring(xml, pretty_print=False, xml_declaration=True, encoding="utf-8",  doctype='<!DOCTYPE TeraStitcher SYSTEM "TeraStitcher.DTD">');
    else:
      return xml;
  else:
    xmlstring = etree.tostring(xml, pretty_print=False, xml_declaration=True, encoding="utf-8",  doctype='<!DOCTYPE TeraStitcher SYSTEM "TeraStitcher.DTD">');
    f = open(xmlFile, 'w')
    f.write(xmlstring);
    f.close();


def xmlImportFile(regularExpression, overlap = None, origin = None, resolution = None, tiling = None, tileExpression = None, zRange = None, xmlFileName = all):
  """Creates the xml import file for TeraStitcher from a directory of image files and optional additional information

  Arguments:
    regularExpression (string): base directory or regular expression to infer the tiled image data, in the latter case use group names
                                names 'col' and 'row' for the tile indices and 'z' for the z-plane indices
    overlap (tuple or None): the (x,y) overlap of the tiles in pixel, if None try to find info from metadata
    origin (tuple or None): the origin of the image, if None then the tpypical value (0.0,0.0,0.0).
    resolution (tupl0 or None): the (x,y,z) resolution of the image in microns per pixel, if None try to find info from metadata
    tiling (tuple or  None): the (row,col) tiling dimensions, if None infer from 'row' and 'col' groups in regularExpression
    tileExpression (function or None): a function of (col,row) that returns a regular expression for the images in the (col,row) tile or None if that tile is missing, 
                                       if None use the regular expression row and col group info
    zRange (function or None. all): a function of (col,row) that returns a the z range specifications for the (col,row) tile
                                    if None use the regular expression z group to infer this from first tile, if all infer this in detail for all tiles             
    xmlFileName (string or None): filename for the xml import file, if None return the xml tree, if all create 'TreaStitcher_import.xml' file in base directory
    
    
  Returns:
    string: filename of the xml import file
    
  Note:
    Example for a regular expression: r'/path/to/image/image_(?P<row>\d{2})_(?P<col>\d{2})_(?P<z>\d{4}).tif'
    
  See also:
    :func:`importData`
  """
  
  ## xml file name
  if xmlFileName is None:
    fh, xmlFileName = tempfile.mkstemp();
    
  ## origin
  if origin is None:
    origin = (0,0,0);
    
    
  #infer size, resolution and overlap
  if size is None or resolution is None or overlap is None:
    firstFile = findFirstFile(regularExpression, sort = False);
    finfo = findFileInfo(firstFile);
    
    if size is None:
      size = finfo['size'];
  
    if resolution is None:
      resolution = finfo['resolution'];
      
    if overlap is None:
      overlap = finfo['overlap'];
    
    for val, name in zip([size, overlap, resolution],['size', 'overlap', 'resolution']) :
      if val is None:
        raise RuntimeError('cannot determine %s from file or input!' % name);
  

  if tiling is None or zRange is None or zRange is all:
    #infer tiling from regular expression
    
    #find file:
    fns, ids = findFileList(regularExpression, groups = ('row', 'col', 'z'));
    ids = np.array(ids);
    
    # get rid of invalid labels
    b = np.zeros(ids.shape[0], dtype = bool);
    for i in range(ids.shape[1]):
      b = np.logical_or(b, np.equal(ids[:,i], None))
 
    fns = fns[~b];
    ids = ids[~b];
    
    if len(fns) == 0:
      raise RuntimeError('no files found that match the expression %s with row, col and z groups' % regularExpression);
    
    # calculate tile dimensions
    rows = np.unique(ids[:,0]);
    nrows = len(rows);
    
    cols = np.unique(ids[:,1]);
    ncols = len(cols);
    
    if tiling is not None and tiling != (nrows, ncols):
      raise RuntimeWarning('specified tiling is different from inferred tiling, min tile number will be used !');
      tiling = (min(nrows, tiling[0]), min(ncols, tiling[1]));
      rows = rows[:tiling[0]];
      cols = cols[:tiling[1]];
    else:
      tiling = (nrows, ncols);
    
    
    # zRanges
    zs  = np.unique(ids[:,2]);
    nzs = len(zs);
    if zRange is None:
      zRange = (0, nzs);
    elif zRange is all:
      zRange = lambda row,col: np.sum(np.logical_and(ids[:,0] == row, ids[:,1] == col));


    #base directory and tile directories
    fnsep = fns[0].split(os.path.sep);
    ffsep = 
    
    #tileExpression
    if tileExpression is None:
      


##############################################################################
### TeraStitcher Runs
##############################################################################

def importData():
  """Runs the import commmand of TeraSticher generating the xml import file
  
  Arguments:
    
  
  
  Returns:
    str: xml import file name
    
  See also:
    :func:`xmlImportFile`
  """
  pass

def alignData():
  """Runs the alignment commmand of TeraSticher aligning the tiles
  
  Arguments:
    
  Returns:
  """
  pass


def stitchData():
  """Runs the alignment commmand of TeraSticher aligning the tiles
  
  Arguments:
    
  Returns:
  """
  pass


def alignData(fixedImage, movingImage, affineParameterFile, bSplineParameterFile = None, resultDirectory = None):
    """Align images using elastix, estimates a transformation :math:`T:` fixed image :math:`\\rightarrow` moving image.
    
    Arguments:
        fixedImage (str): image source of the fixed image (typically the reference image)
        movingImage (str): image source of the moving image (typically the image to be registered)
        affineParameterFile (str or None): elastix parameter file for the primary affine transformation
        bSplineParameterFile (str or None): elastix parameter file for the secondary non-linear transformation
        resultDirectory (str or None): elastic result directory
        
    Returns:
        str: path to elastix result directory
    """
    
    checkElastixInitialized();
    global ElastixBinary;
    
    if resultDirectory == None:
        resultDirectory = tempfile.gettempdir();
    
    if not os.path.exists(resultDirectory):
        os.mkdir(resultDirectory);
    
    
    if bSplineParameterFile is None:
        cmd = ElastixBinary + ' -threads 16 -m ' + movingImage + ' -f ' + fixedImage + ' -p ' + affineParameterFile + ' -out ' + resultDirectory;
    elif affineParameterFile is None:
        cmd = ElastixBinary + ' -threads 16 -m ' + movingImage + ' -f ' + fixedImage + ' -p ' + bSplineParameterFile + ' -out ' + resultDirectory;
    else:
        cmd = ElastixBinary + ' -threads 16 -m ' + movingImage + ' -f ' + fixedImage + ' -p ' + affineParameterFile + ' -p ' + bSplineParameterFile + ' -out ' + resultDirectory;
        #$ELASTIX -threads 16 -m $MOVINGIMAGE -f $FIXEDIMAGE -fMask $FIXEDIMAGE_MASK -p  $AFFINEPARFILE -p $BSPLINEPARFILE -out $ELASTIX_OUTPUT_DIR
    
    res = os.system(cmd);
    
    if res != 0:
        raise RuntimeError('alignData: failed executing: ' + cmd);
    
    return resultDirectory


def transformData(source, sink = [], transformParameterFile = None, transformDirectory = None, resultDirectory = None):
    """Transform a raw data set to reference using the elastix alignment results
    
    If the map determined by elastix is 
    :math:`T \\mathrm{fixed} \\rightarrow \\mathrm{moving}`, 
    transformix on data works as :math:`T^{-1}(\\mathrm{data})`.
        
    Arguments:
        source (str or array): image source to be transformed
        sink (str, [] or None): image sink to save transformed image to. if [] return the default name of the data file generated by transformix.
        transformParameterFile (str or None): parameter file for the primary transformation, if None, the file is determined from the transformDirectory.
        transformDirectory (str or None): result directory of elastix alignment, if None the transformParameterFile has to be given.
        resultDirectory (str or None): the directorty for the transformix results
        
    Returns:
        array or str: array or file name of the transformed data
    """
    
    global TransformixBinary;    
    
    if isinstance(source, numpy.ndarray):
        imgname = os.path.join(tempfile.gettempdir(), 'elastix_input.tif');
        io.writeData(source, imgname);
    elif isinstance(source, basestring):
        if io.dataFileNameToType(source) == "TIF":
            imgname = source;
        else:
            imgname = os.path.join(tempfile.gettempdir(), 'elastix_input.tif');
            io.transformData(source, imgname);
    else:
        raise RuntimeError('transformData: source not a string or array');

    if resultDirectory == None:
        resultdirname = os.path.join(tempfile.tempdir, 'elastix_output');
    else:
        resultdirname = resultDirectory;
        
    if not os.path.exists(resultdirname):
        os.makedirs(resultdirname);
        
    
    if transformParameterFile == None:
        if transformDirectory == None:
            raise RuntimeError('neither alignment directory and transformation parameter file specified!'); 
        transformparameterdir = transformDirectory
        transformParameterFile = getTransformParameterFile(transformparameterdir);
    else:
        transformparameterdir = os.path.split(transformParameterFile);
        transformparameterdir = transformparameterdir[0];
    
    #transform
    #make path in parameterfiles absolute
    setPathTransformParameterFiles(transformparameterdir);
   
    #transformix -in inputImage.ext -out outputDirectory -tp TransformParameters.txt
    cmd = TransformixBinary + ' -in ' + imgname + ' -out ' + resultdirname + ' -tp ' + transformParameterFile;
    
    res = os.system(cmd);
    
    if res != 0:
        raise RuntimeError('transformData: failed executing: ' + cmd);
    
    
    if not isinstance(source, basestring):
        os.remove(imgname);

    if sink == []:
        return getResultDataFile(resultdirname);
    elif sink is None:
        resultfile = getResultDataFile(resultdirname);
        return io.readData(resultfile);
    elif isinstance(sink, basestring):
        resultfile = getResultDataFile(resultdirname);
        return io.convertData(resultfile, sink);
    else:
        raise RuntimeError('transformData: sink not valid!');


 
##############################################################################
### Test
##############################################################################

     
def test():
    """Test Stiching module"""
    import ClearMap.Alignment.Elastix as self
    reload(self)
    
    from ClearMap.Settings import ClearMapPath;
    import os, numpy
    
    p = ClearMapPath;
    
    resultdir = os.path.join(p, 'Test/Elastix/Output');
    
    print 'Searching for transformation parameter file in ' + resultdir;
    pf = self.getTransformParameterFile(resultdir)
      
    print 'Found: ' + pf;
    
    
    #replace path in trasform parameter files:
    self.setPathTransformParameterFiles(resultdir)
    
    #initialize
    self.initializeElastix('/home/ckirst/programs/elastix')
    self.printSettings()

    #transform points
    pts = numpy.random.rand(5,3);    
     
    print 'Transforming points: '
    tpts = self.transformPoints(pts, transformParameterFile = pf, indices = False);
    print pts
    print 'Transformed points: '
    print tpts
    
    
    #deformation and distance fields     
    df = self.deformationField(transformParameterFile = pf, resultDirectory = None);
    #df = '/tmp/elastix_output/deformationField.mhd';

    import ClearMap.IO as io
    data = io.readData('/tmp/elastix_output/deformationField.mhd');
    
    ds = self.deformationDistance(data);
    
    io.writeData(os.path.join(p, 'Test/Elastix/Output/distances.raw'), ds);
    
if __name__ == "__main__":
    test();
    





  

    
    
    

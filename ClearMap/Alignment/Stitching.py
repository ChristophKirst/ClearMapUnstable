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
import glob
import shutil
import tempfile
import numpy as np
import re
import natsort


from lxml import etree


import ClearMap.Settings as settings
import ClearMap.IO as io
import ClearMap.Visualization.Plot as cplt;

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
    
    return np.array(files), np.array(infos);



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
    search = re.compile('^' + n + '$').search;
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
    imgxml = etree.fromstring(str(imgxml));
    
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


def baseDirectoryFromXMLImport(xmlFile):
  """Extracts the base directory specified in the xml descriptor
  
  Arguments:
    xmlFile: the TeraStitcher xml file
  
  Returns:
    str: base directory specified int he xml file
  """
  
  try:
    xmlt = etree.parse(xmlFile);
    baseDirectory = xmlt.find('stacks_dir').get('value');
    return baseDirectory;
  except:
    return None;
 

def clearBaseDirectory(baseDirectory):
  """Clears the base directory from TeraStitcher generated binary and tif files
  
  Arguments:
    baseDirectory (str): the base directory or xml descriptor file
    
  Note:
    TeraStitcher produces a 'mdata.bin' and 'test_middle_slice.tif' file during import
  """
  
  bD = baseDirectoryFromXMLImport(baseDirectory);
  if bD is None:
    bD = baseDirectory;

  try:
    binfile = os.path.join(baseDirectory, 'mdata.bin');
    print 
    os.remove(binfile);
  except:
    pass;
  
  try:
    testfile = os.path.join(baseDirectory, 'test_middle_slice.tif');
    os.remove(testfile);
  except:
    pass;
  
  
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
  
  s.set('ABS_V', str(int(displacement[1])));
  s.set('ABS_H', str(int(displacement[0])));
  if len(displacement) > 2:
    s.set('ABS_D', str(int(displacement[2]))); 
  else:
    s.set('ABS_D', str(0));
  
  if stitchable:
    s.set('STITCHABLE', 'yes');
  else:
    s.set('STITCHABLE', 'no');
  
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


def xmlImport(baseDirectory, size, resolution = (1.0, 1.0, 1.0), origin = (0.0, 0.0, 0.0), overlap = (0.0, 0.0), 
              tiling = (1,1), tileExpression = None, zRange = None, rows = None, cols = None, form = None,
              xmlImportFile = None, asString = False):
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
                                       images in the tile directory are used
    rows,cols (list or None): if list use these values for row and col indices, if None integer ranges are used
    zRange (function or None): zRange(row, col) should return the zrange for 
                               the (row,col) tile as a string, if None, full 
                               range is used
    form (str or None): volume format, if none use 'TiledXY|2Dseries'
    xmlImportFile (str or None): file name to save the xml to, if None return xml 
                               string or xml element tree
    asString (bool): if True return as string otherwise as xml element tree
  """
  
  xml = etree.Element('TeraStitcher');
  
  if form is None:
    form = 'TiledXY|2Dseries';
  xml.set('volume_format', form);
  
  e = etree.Element('stacks_dir');
  e.set('value', baseDirectory);  
  xml.append(e);
  
  e = etree.Element('voxel_dims');
  e.set('V', str(resolution[1]));
  e.set('H', str(resolution[0]));
  e.set('D', str(resolution[2]));
  xml.append(e);
  
  e = etree.Element('origin');
  e.set('V', str(origin[1]));
  e.set('H', str(origin[0]));
  e.set('D', str(origin[2]));
  xml.append(e);

  e = etree.Element('mechanical_displacements');
  displ = displacements(size, resolution, overlap, units = 'Microns');
  e.set('V', str(displ[1]));
  e.set('H', str(displ[0]));
  xml.append(e);
  
  e = etree.Element('dimensions');
  e.set('stack_columns', str(tiling[0]));
  e.set('stack_rows', str(tiling[1]));
  e.set('stack_slices', str(size[2]));
  xml.append(e);
  
  if zRange is all or zRange is None:
    zRange = lambda row,col: (0, size[2]); 
  if isinstance(zRange, tuple) or isinstance(zRange, list):
    zr = zRange;
    zRange = lambda row,col: zr
  
  displ = displacements(size, resolution, overlap, units = 'Pixel');
    
  if tileExpression is None:
    tileExpression = lambda row,col: '';
  
  stack = etree.Element('STACKS');
  
  if rows is None:
    rows = range(tiling[0]);
  if cols is None:
    cols = range(tiling[1]);
  
  for rowi, row in enumerate(rows):
    for coli, col in enumerate(cols):
      stackdispl = [displ[0] * coli, displ[1] * rowi];
      te = tileExpression(row,col);
      dn,fn = os.path.split(te);      
      e = xmlImportStack(rowi, coli, zRange(row,col), stackdispl, stitchable = False, directory = dn, expression = fn);
      stack.append(e);
  
  xml.append(stack);
  
  
  if xmlImportFile is None:
    if asString:
      return etree.tostring(xml, pretty_print=True, xml_declaration=True, encoding="utf-8",  doctype='<!DOCTYPE TeraStitcher SYSTEM "TeraStitcher.DTD">');
    else:
      return xml;
  else:
    xmlstring = etree.tostring(xml, pretty_print=True, xml_declaration=True, encoding="utf-8",  doctype='<!DOCTYPE TeraStitcher SYSTEM "TeraStitcher.DTD">');
    f = open(xmlImportFile, 'w')
    f.write(xmlstring);
    f.close();
    return xmlImportFile;


def xmlImportFile(regularExpression, size = None, overlap = None, origin = None, resolution = None, tiling = None, tileExpression = None, zRange = None, xmlImportFile = None, asString = False):
  """Creates the xml import file for TeraStitcher from a directory of image files and optional additional information

  Arguments:
    regularExpression (string): base directory or regular expression to infer the tiled image data, in the latter case use group names
                                names 'col' and 'row' for the tile indices and 'z' for the z-plane indices
    size (tuple or None): volume size in pixel if None inferred from meta data
    overlap (tuple or None): the (x,y) overlap of the tiles in pixel, if None try to find info from metadata
    origin (tuple or None): the origin of the image, if None then the tpypical value (0.0,0.0,0.0).
    resolution (tupl0 or None): the (x,y,z) resolution of the image in microns per pixel, if None try to find info from metadata
    tiling (tuple or  None): the (row,col) tiling dimensions, if None infer from 'row' and 'col' groups in regularExpression
    tileExpression (function or None): a function of (col,row) that returns a regular expression for the images in the (col,row) tile or None if that tile is missing, 
                                       if None use the regular expression row and col group info
    zRange (function or None. all): a function of (col,row) that returns a the z range specifications for the (col,row) tile
                                    if None use the regular expression z group to infer this from first tile, if all infer this in detail for all tiles             
    xmlImportFile (string or None): filename for the xml import file, if None return the xml tree, if all create 'TreaStitcher_import.xml' file in base directory
    asString (bool): if True return xml code as string, otherwise return as a xml tree
    
  Returns:
    string: filename of the xml import file
    
  Note:
    Example for a regular expression: r'/path/to/image/image_(?P<row>\d{2})_(?P<col>\d{2})_(?P<z>\d{4}).tif'
    
  See also:
    :func:`importData`
  """
  
  ## origin
  if origin is None:
    origin = (0,0,0);
  
  #infer size, resolution and overlap
  if size is None or resolution is None or overlap is None:
    firstFile = findFirstFile(regularExpression);
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
    b = np.logical_not(b);
    fns = fns[b];
    ids = ids[b];
    
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
      zRange = lambda row,col: (0, np.sum(np.logical_and(ids[:,0] == row, ids[:,1] == col)));
    
  else:
    rows = None;
    cols = None;
    fns = None;
    
    nzs = 0;
    for row in range(tiling[0]):
      for col in range(tiling[1]):
        nzs = max(nzs, zRange(row,col)[1]);
    
  size = tuple(size) + (nzs,);
  
  #base directory and tile directories
  if fns is None:
    if firstFile is None:
      fn = findFirstFile(regularExpression, sort = False);
    else:
      fn = firstFile;
  else:
    fn = fns[0];
  
  fnsep = fn.split(os.path.sep);
  fesep = regularExpression.split(os.path.sep);
  
  if len(fnsep) != len(fesep):
    raise RuntimeError('inconsistent file names and file expression!');
  
  for i in range(len(fnsep)):
    if fnsep[i] != fesep[i]:
      baseDirectory = os.path.sep.join(fesep[:i]);
      regularExpression = os.path.sep.join(fesep[i:]);
      break;
  
  #tileExpression
  if tileExpression is None:
    def makeTileExpression(row,col):
      te = re.sub(r'\(\?\P\<row\>.*?\)', str(row), regularExpression, count = 1);
      return re.sub(r'\(\?\P\<col\>.*?\)', str(col), te, count = 1);
    tileExpression = makeTileExpression;
  
  # create xml import 
  return xmlImport(baseDirectory, size = size, resolution = resolution, origin = origin, overlap = overlap, 
            tiling = tiling, tileExpression = tileExpression, zRange = zRange, rows = rows, cols = cols,
            xmlImportFile = xmlImportFile, asString = asString);
    

##############################################################################
### TeraStitcher Runs
##############################################################################

def importData(xmlImportFile = None, baseDirectory = None, resolution = None, orientation = None,
               regularExpression = None, form = None, rescan = None, sparse = None, xmlResultFile = None, clear = True):
  """Runs the import commmand of TeraSticher generating the xml import file from a folder of files
  
  Arguments:
    xmlImportFile (str or None): the xml import descriptor
    baseDirectory (str or None): the base directory of the image data, None if xmlImportFile is given
    regularExpression (str or None): optional regular expression for images
    resolution (tuple or None): optional resolution of the image in micorns per pixel
    orientation (tuple or None): optional orientation of the image 
    form (str or None): the import format (if None the default 'TiledXY|2Dseries' is used)
    rescan (bool or None): optional rescan of image data
    sparse (bool or None): optional switch to specify if data is sparse
    xmlResultFileName (str or None): output xml file name
    clear (bool): if True delete previous imported binary files
  
  Returns:
    str: xml import file name
    
  See also:
    :func:`xmlImportFile`, `TeraStitcher import step <https://github.com/abria/TeraStitcher/wiki/Step-1:-Import>`_.
  """
  
  checkSticherInitialized();
  global TeraStitcherBinary;
    
  if baseDirectory is None and xmlImportFile is None:
    raise RuntimeError('importData requires baseDirectory or xmlImportFile!');
  
  if xmlImportFile is None:
    if resolution is None:
      resolution = (1.0, 1.0, 1.0);
    if orientation is None:
      orientation = (1,2,3);
  
  cmd = TeraStitcherBinary + ' --import ';
  
  if baseDirectory is not None:
    cmd = cmd + ' --volin="' + baseDirectory + '" ';
  
  if xmlImportFile is not None:
    cmd = cmd + ' --projin="' + xmlImportFile + '" ';
  

  #voxel size / resolution
  if resolution is not None:
    vsize = ['--vxl1=', '--vxl2=', '--vxl3='];  
    for i in range(3):
      cmd = cmd + vsize[i] + str(resolution[i]) + ' ';
  
  #reference orientation1
  if orientation is not None:
    ref = ['--ref1=', '--ref2=', '--ref3='];  
    #refnames = ['x', 'y', 'z'];
    for i in range(3):
      #if orientation[i] < 0:
      #  rn = '-'
      #else:
      #  rn = '';
      #rn = rn + refnames[orientation[abs(i)]];
      rn = str(orientation[i]);
      cmd = cmd + ref[i] + rn + ' ';
  
  #optional arguments
  if xmlResultFile is not None:
    cmd = cmd + '--projout="' + xmlResultFile + '" ';
  else:
    if xmlImportFile is not None:
      xmlResultFile = xmlImportFile;
    else:
      xmlResultFile = 'xml_import.xml';
  
  if form is not None:
    cmd = cmd + '--volin_plugin="' + form + '" ';
  #else:
  #  cmd = cmd + '--volin_plugin="TiledXY|2Dseries" ';
  
  if rescan is True:
    cmd = cmd + '--rescan ';
    
  if sparse is True:
    cmd = cmd + '--sparse_data ';
    
  if regularExpression is not None:
    cmd = cmd + '--imin_regex="' + regularExpression + '" ';
  
  
  if clear:
    if baseDirectory is None:
      baseDirectory = baseDirectoryFromXMLImport(xmlImportFile);
    clearBaseDirectory(baseDirectory);
  
  print 'running: %s' % cmd
  res = os.system(cmd);
  
  if res != 0:
    raise RuntimeError('importData: failed executing: ' + cmd);
    
  return xmlResultFile;


def readImportPreview(importResult):
  """Returns the prestiched preview file to check orientations and coarse alignment
  
  Arguments:
    importResult (str): the base directory of the image data or the import xml file
  
  Returns:
    array: preview image
  """
  
  baseDirectory = baseDirectoryFromXMLImport(importResult) 
  if baseDirectory is None:
    baseDirectory = importResult;

  fn = os.path.join(baseDirectory, 'test_middle_slice.tif');
  return io.readData(fn);
  

def plotImportPreview(importResult):
  """Plots the prestiched preview file to check orientations and coarse alignment
  
  Arguments:
    importResult (str): the base directory of the image data or the import xml file
  
  Returns:
    object: plot axes
  """
  return cplt.plotTiling(readImportPreview(importResult));

  
def alignData(xmlImportFile, slices = None, subRegion = None, overlap = None, search = None, algorithm = None, channel = None, xmlResultFile = None):
  """Runs the alignment commmand of TeraSticher aligning the tiles
  
  Arguments:
    xmlImportFile (str or None): the xml import descriptor
    slices (int or None): optional number of top slices to use for alignment
    subRegion (tuple or None): optional sub region in the form ((xmin,xmax),(ymin,ymax), (zmin, zmax))
    overlap (tuple or None): optional estimated overlap
    search (tuple or None): optional search region as pixel radii, i.e. (xradius, yradius)
    algorithm (str or None): optional algorithm 'MIPNCC' for MIP Normalized Cross-Correlation or 'PC' for Phase Correlation 
    channel (str or None): optional channel to use 'R', 'G', 'B' or 'all'
    xmlResultFile (str or None): output xml displacement file name
  
  Returns:
    str : the result xml displacement file name
    
  See also:
    `TeraStitcher align step <https://github.com/abria/TeraStitcher/wiki/Step-2:-Align>`_.
  """
  
  checkSticherInitialized();
  global TeraStitcherBinary;
  
  cmd = TeraStitcherBinary + ' --displcompute ';
  
  cmd = cmd + ' --projin="' + xmlImportFile + '" ';
  
  #slices
  if slices is not None:
    cmd = cmd + ' --subvoldim=' + str(slices) + ' ';
  
  if subRegion is not None:
    sns = (('--R0=', '--R1='), ('--C0=', '--C1='), ('--D0=', '--D1-'));
    for d in range(3):
      for m in range(2):
        if subRegion[d][m] is not None:
          cmd = cmd + sns[d][m] + str(subRegion[d][m]) + ' ';
    
  if overlap is not None:
    ov = ('--oV=', '--oH=');
    for d in range(2):
      if overlap[d] is not None:
        cmd = cmd + ov[d] + str(overlap[d]) + ' ';
        
  if search is not None:
    sd = ('--sV=', '--sH=', '--sD=');
    for d in range(3):
      if search[d] is not None:
        cmd = cmd + sd[d] + str(search[d]) + ' ';
    
  if algorithm is not None:
     cmd = cmd + '--algorithm="' + algorithm + '" ';
  
  if channel is not None:
    cmd = cmd + '--imin_channel="' + channel + '" ';
  
  if xmlResultFile is not None:
    cmd = cmd + '--projout="' + xmlResultFile + '" ';
  else:
    xmlResultFile = 'xml_displcomp.xml';
  
  print 'running ' + cmd
  res = os.system(cmd);
  
  if res != 0:
    raise RuntimeError('alignData: failed executing: ' + cmd);
    
  return xmlResultFile;


def projectDisplacements(xmlDisplacementFile, xmlResultFile = None):
  """Runs the projection step of TeraSticher to aligning the tiles
  
  Arguments:
    xmlDisplacementFile (str or None): the xml displacement descriptor
    xmlResultFile (str or None): output xml projection file name
  
  Returns:
    str : the result xml projection file name
    
  See also:
    `TeraStitcher project step <https://github.com/abria/TeraStitcher/wiki/Step-3:-Project>`_.
  """
  
  checkSticherInitialized();
  global TeraStitcherBinary;
  
  cmd = TeraStitcherBinary + ' --displproj ';
  
  cmd = cmd + ' --projin="' + xmlDisplacementFile + '" ';
  
  if xmlResultFile is not None:
    cmd = cmd + '--projout="' + xmlResultFile + '" ';
  else:
    xmlResultFile = 'xml_displproj.xml';
  
  print 'running: ' + cmd;
  res = os.system(cmd);
  
  if res != 0:
    raise RuntimeError('projectDisplacements: failed executing: ' + cmd);
    
  return xmlResultFile;


def thresholdDisplacements(xmlProjectionFile, threshold = None, xmlResultFile = None):
  """Runs the thresholding step of TeraSticher to aligning the tiles
  
  Arguments:
    xmlProjectionFile (str or None): the xml projection descriptor
    threshold (float or None): optional threshold value
    xmlResultFile (str or None): output xml thresholded file name
  
  Returns:
    str : the result xml thresholded file name
    
  See also:
    `TeraStitcher project step <https://github.com/abria/TeraStitcher/wiki/Step-4:-Threshold>`_.
  """
  
  checkSticherInitialized();
  global TeraStitcherBinary;
  
  cmd = TeraStitcherBinary + ' --displthres ';
  
  cmd = cmd + ' --projin="' + xmlProjectionFile + '" ';
  
  if threshold is None:
    threshold = 0;
  cmd = cmd + ' --threshold=' + str(threshold) + ' ';
  
  if xmlResultFile is not None:
    cmd = cmd + '--projout="' + xmlResultFile + '" ';
  else:
    xmlResultFile = 'xml_displthres.xml';
  
  print 'running: ' + cmd;
  res = os.system(cmd);
  
  if res != 0:
    raise RuntimeError('thresholdDisplacements: failed executing: ' + cmd);
    
  return xmlResultFile;




def placeTiles(xmlThresholdFile, algorithm = None, xmlResultFile = None):
  """Runs the placement step of TeraSticher to aligning the tiles
  
  Arguments:
    xmlThresholdFile (str or None): the xml threshold descriptor
    algorithm (str or None): optional algorithm to use for placement: 
                             'MST' for Minimum Spanning Tree
                             'SCAN_V' for scanning along vertical axis
                             'SCAN_H' for scannning along horizontal axis
    xmlResultFile (str or None): output xml placed file name
  
  Returns:
    str : the result xml placed file name
    
  See also:
    `TeraStitcher project step <https://github.com/abria/TeraStitcher/wiki/Step-5:-Place>`_.
  """
  
  checkSticherInitialized();
  global TeraStitcherBinary;
  
  cmd = TeraStitcherBinary + ' --placetiles ';
  
  cmd = cmd + ' --projin="' + xmlThresholdFile + '" ';
  
  if algorithm is not None:
    cmd = cmd + ' --algorithm="' + algorithm + '" ';
  
  if xmlResultFile is not None:
    cmd = cmd + '--projout="' + xmlResultFile + '" ';
  else:
    xmlResultFile = 'xml_merging.xml';
  
  print 'running: ' + cmd;
  res = os.system(cmd);
  
  if res != 0:
    raise RuntimeError('placeTiles: failed executing: ' + cmd);
    
  return xmlResultFile;


def stitchData(xmlPlacementFile, resultPath, algorithm = None, resolutions = None, form = None, channel = None, subRegion = None, bitDepth = None, blockSize = None, clear = True):
  """Runs the final stiching step of TeraSticher
  
  Arguments:
    xmlPlacementFile (str or None): the xml placement descriptor
    resultPath (str): result path or file name for the stiched data
    algorithm (str or None): optional algorithm to use for placement: 
                             'NO_BLEND' for no blending
                             'SIN_BLEND' for sinusoidal blending
                             'STACK_MARGIN' for show stack margin
    resolutions (tuple or None): the different resolutions to produce
    form (str or None): the output form, if None determined automatically
    channel (str or None): the channels to use, 'R', 'G', 'B' or 'all'
    subRegion (tuple or None): optional sub region in the form ((xmin,xmax),(ymin,ymax), (zmin, zmax))
    bitDepth (int or None): the pits per pixel to use, default is 8
    blockSize (tuple): the sizes of various blocks to save stiched image into
  
  Returns:
    str : the result path or file name of the stiched data
    
  See also:
    `TeraStitcher project step <https://github.com/abria/TeraStitcher/wiki/Step-6:-Merge>`_.
  """
  
  checkSticherInitialized();
  global TeraStitcherBinary;
  
  cmd = TeraStitcherBinary + ' --merge --imout_format="tif" ';
  
  cmd = cmd + '--projin="' + xmlPlacementFile + '" ';
  
  if len(resultPath) > 3 and resultPath[-4:] == '.tif':
    form = 'TiledXY|3Dseries';
    resultPath, filename = os.path.split(resultPath);
  else:
    filename = None;
  
  cmd = cmd + '--volout="' + resultPath + '" ';
  
  if algorithm is not None:
    cmd = cmd + ' --algorithm="' + algorithm + '" ';
    
  if resolutions is not None:
    cmd = cmd + '--resolutions="'
    for r in sorted(resolutions):
      cmd = cmd + str(r);
    cmd = cmd + ' ';
    
  if form is not None:
    cmd = cmd + '--volout_plugin="' + form + '" ';
  
  if channel is not None:
    cmd = cmd + '--imin_channel="' + channel + '" ';
  
  if subRegion is not None:
    sns = (('--R0=', '--R1='), ('--C0=', '--C1='), ('--D0=', '--D1-'));
    for d in range(3):
      for m in range(2):
        if subRegion[d][m] is not None:
          cmd = cmd + sns[d][m] + str(subRegion[d][m]) + ' ';
          
  if blockSize is not None:
    bs = ('--slicewidth=', '--sliceheight=', '--slicedepth=');
    for d in range(3):
      if blockSize[d] is not None:
        cmd = cmd + bs[d] + str(blockSize[d]) + ' ';
  
  if bitDepth is not None:
    cmd = cmd + '--imout_depth=' + str(bitDepth) + ' ';
 
  print 'running: ' + cmd;
  res = os.system(cmd);
  
  if res != 0:
    raise RuntimeError('stitchData: failed executing: ' + cmd);
  
  if filename is not None:
    #get most recent created file and move to desired filename
    imgfile = max(glob.glob(os.path.join(resultPath, '*/*/*/*')), key = os.path.getmtime);
    filename = os.path.join(resultPath, filename);
    os.rename(imgfile, filename);
    if clear:
      imgpath = os.path.sep.join(imgfile.split(os.path.sep)[:-3]);
      shutil.rmtree(imgpath)
    return filename;
  else:
    return resultPath;

##############################################################################
### Test
##############################################################################

     
def test():
    """Test Stiching module"""
    import ClearMap.Alignment.Stitching as st
    reload(st)

    st.printSettings()
    
    from ClearMap.Settings import ClearMapPath;
    import os
    
    datadir = os.path.join(ClearMapPath, 'Test/Data/Stitching');
    outdir = os.path.join(ClearMapPath, 'Test/Data/Stitching/Output');
    
    #import glob
    #fn = os.path.join(datadir, '14-22-21_TH-s3-ov45-na004_UltraII[[]* x *[]]_C00_xyz-Table Z*.ome.tif')
    #fn = os.path.join(datadir, '*')
    #print fn
    #print glob.glob(fn)
    
    fileExpression = os.path.join(datadir, r'14-22-21_TH-s3-ov45-na004_UltraII\[(?P<row>\d{2}) x (?P<col>\d{2})\]_C00_xyz-Table Z(?P<z>\d{4}).ome.tif')
    
    st.findFileList(fileExpression)
    
    ff = st.findFirstFile(fileExpression, sort = False)
    
    fi = st.findFileInfo(ff);
    print fi
    
    ### import file
    reload(st); 
    #xmlstring = st.xmlImportFile(fileExpression, overlap = None, origin = None, resolution = None, tiling = None, tileExpression = None, zRange = None, xmlImportFile = None, asString = True)
    #print xmlstring
    
    importfile = st.xmlImportFile(fileExpression, xmlImportFile = os.path.join(datadir, 'TeraStitcher_import.xml'));
    print importfile

    ### run import
    st.importData(importfile)
    
    ### preview
    st.plotImportPreview(importfile)
    
    ### run alignment chain
    reload(st)
    alignfile = st.alignData(importfile, slices = 5, search = (10,10,0), xmlResultFile = os.path.join(datadir, 'TeraStitcher_align.xml'));

    ### project
    reload(st);
    projectfile = st.projectDisplacements(alignfile,  xmlResultFile = os.path.join(datadir, 'TeraStitcher_project.xml'))
    
    
    ### threshold
    reload(st);
    thresfile = st.thresholdDisplacements(projectfile,  xmlResultFile = os.path.join(datadir, 'TeraStitcher_threshold.xml'))

    ### place
    reload(st);
    placefile = st.placeTiles(thresfile, xmlResultFile = os.path.join(datadir, 'TeraStitcher_place.xml'))
    
    ### stitch it
    reload(st);
    result = st.stitchData(placefile, resultPath = os.path.join(datadir, 'stiched.tif'))
    result = st.stitchData(placefile, resultPath = os.path.join(datadir, 'stiched.tif'))
    
if __name__ == "__main__":
    test();
    





  

    
    
    

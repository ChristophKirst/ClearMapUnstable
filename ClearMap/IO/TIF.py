# -*- coding: utf-8 -*-
"""
Interface to tif image files and stacks.

The interface makes use of the tifffile library.

Example:
    >>> import os, numpy
    >>> import ClearMap.IO.TIF as tif
    >>> from ClearMap.Settings import ClearMapPath
    >>> filename = os.path.join(ClearMapPath,'Test/Data/Tif/composite.tif') 
    >>> data = tif.readData(filename);
    >>> print data.shape
    (50, 100, 30, 4)
    
"""
#:copyright: Copyright 2015 by Christoph Kirst, The Rockefeller University, New York City
#:license: GNU, see LICENSE.txt for details.



import numpy
import tifffile as tiff

import ClearMap.IO as io


def dataSize(filename, **args):
    """Returns size of data in tif file
    
    Arguments:
        filename (str): file name as regular expression
        x,y,z (tuple): data range specifications
    
    Returns:
        tuple: data size
    """
    
    t = tiff.TiffFile(filename);
    d3 = len(t.pages);
    d2 = t.pages[0].shape;
    #d2 = (d2[0], d2[1]);
    d2 = (d2[1], d2[0]);
    if d3 > 1:
        dims = d2 + (d3,);
    else:
        dims =  d2;
    
    return io.dataSizeFromDataRange(dims, **args);

def dataZSize(filename, z = all, **args):
    """Returns z size of data in tif file
    
    Arguments:
        filename (str): file name as regular expression
        z (tuple): z data range specification
    
    Returns:
        int: z data size
    """
    
    t = tiff.TiffFile(filename);
    d3 = len(t.pages);
    if d3 > 1:
        return io.toDataSize(d3, r = z);
    else:
        return None;



def readData(filename, x = all, y = all, z = all, **args):
    """Read data from a single tif image or stack
    
    Arguments:
        filename (str): file name as regular expression
        x,y,z (tuple): data range specifications
    
    Returns:
        array: image data
    """
    
    dsize = dataSize(filename);
    #print "dsize %s" % str(dsize);    
    
    if len(dsize) == 2:
        data = tiff.imread(filename, key = 0);
        #print "data.shape %s" % str(data.shape);        
        
        return io.dataToRange(data.transpose([1,0]), x = x, y = y);
        #return io.dataToRange(data, x = x, y = y);
        
    else:
        if z is all:
            data = tiff.imread(filename);
            if data.ndim == 2:
                # data = data
                data = data.transpose([1,0]);
            elif data.ndim == 3:
                #data = data.transpose([1,2,0]);
                data = data.transpose([2,1,0]);
            elif data.ndim == 4: # multi channel image
                #data = data.transpose([1,2,0,3]);
                data = data.transpose([2,1,0,3]);
            else:
                raise RuntimeError('readData: dimension %d not supproted!' % data.ndim)
            
            return io.dataToRange(data, x = x, y = y, z = all);
        
        else: #optimize for z ranges
            ds = io.dataSizeFromDataRange(dsize, x = x, y = y, z = z);
            t = tiff.TiffFile(filename);
            p = t.pages[0];
            data = numpy.zeros(ds, dtype = p.dtype);
            rz = io.toDataRange(dsize[2], r = z);
            
            #print "test"
            #print rz;
            #print dsize            
            
            for i in range(rz[0], rz[1]):
                xydata = t.pages[i].asarray();
                #data[:,:,i-rz[0]] = io.dataToRange(xydata, x = x, y = y);
                data[:,:,i-rz[0]] = io.dataToRange(xydata.transpose([1,0]), x = x, y = y);
            
            return data


def writeData(filename, data):
    """Write image data to tif file
    
    Arguments:
        filename (str): file name 
        data (array): image data
    
    Returns:
        str: tif file name
    """
    
    d = len(data.shape);
    
    if d == 2:
        #tiff.imsave(filename, data);
        tiff.imsave(filename, data.transpose([1,0]));
    elif d == 3:   
        #tiff.imsave(filename, data.transpose([2,0,1]));
        tiff.imsave(filename, data.transpose([2,1,0]));
    elif d == 4:
        #tiffile (z,y,x,c)
        t = tiff.TiffWriter(filename, bigtiff = True);
        #t.save(data.transpose([2,0,1,3]), photometric = 'minisblack',  planarconfig = 'contig');
        t.save(data.transpose([2,1,0,3]), photometric = 'minisblack',  planarconfig = 'contig')
        t.close();    
    else:
        raise RuntimeError('writing multiple channel data to tif not supported');
    
    return filename;
    

def copyData(source, sink):
    """Copy a data file from source to sink
    
    Arguments:
        source (str): file name pattern of source
        sink (str): file name pattern of sink
    
    Returns:
        str: file name of the copy
    """ 
    
    io.copyFile(source, sink);


def readMetaData(source, info = all):
  """Reads the meta data from the image file
  
  Arguments:
    source: the data source
    info (list or all): optional list of keywords
  
  Returns:
    object: an object with the meta data
  """
  
  from PIL import Image
  from PIL.ExifTags import TAGS
  from lxml import etree
     
  def ome_info(fn):
      ret = {}
      i = Image.open(fn)
      info = i.tag.as_dict();
      for tag, value in info.iteritems():
          decoded = TAGS.get(tag)
          ret[decoded] = value
      return ret
      
  imginfo = ome_info(source);
  
  if info is all:
    return imginfo;
    
  
  keys = imginfo.keys();
  subinfo = {x : None for x in info};

  # get image sizes    
  if 'size' in info:
    subinfo['size'] = dataSize(source);
  
  if 'resolution' in info or 'overlap' in info:
    try:
      imgxml = imginfo['ImageDescription'];
      if isinstance(imgxml, tuple):
        imgxml = imgxml[0];
      imgxml = etree.fromstring(str(imgxml));
    except:
      return subinfo;
    
    # get resolution -> ome generic
    if 'resolution' in info:        
      try:
        pix = [x for x in imgxml.iter('{*}Pixels')];
        pix = pix[0].attrib;
        keys = pix.keys();
        if 'PhysicalSizeX' in keys and 'PhysicalSizeY' in keys and 'PhysicalSizeZ' in keys:
          subinfo['resolution'] = (float(pix['PhysicalSizeX']), float(pix['PhysicalSizeY']), float(pix['PhysicalSizeZ']));
      except:
        pass
  
    # get overlap -> not ome generic ???
    if 'overlap' in info:
      try:
        e1 = [x for x in imgxml.iter('{*}xyz-Table_X_Overlap')];
        e2 = [x for x in imgxml.iter('{*}xyz-Table_Y_Overlap')];
        subinfo['overlap'] = (float(e1[0].attrib['Value']), float(e2[0].attrib['Value']));
      except:
        pass
      
  return subinfo;



def test():    
    """Test TIF module"""  
    import ClearMap.IO.TIF as self
    reload(self)
    
    from ClearMap.Settings import ClearMapPath as basedir
    import os
    import numpy

    fn = os.path.join(basedir,'Test/Data/Tif/test.tif')  

    data = numpy.random.rand(20,50,10);
    data[5:15, 20:45, 2:9] = 0;
    data = 20 * data;
    data = data.astype('int32');

    #reload(self)
    print "writing raw image to: " + fn;    
    self.writeData(fn, data);

    print "Loading raw image from: " + fn;
    img = self.readData(fn);  
    print "Image size: " + str(img.shape)
    
    diff = img - data;
    print (diff.max(), diff.min())
    
    
    print "Loading raw image from %s with limited z range: " % fn;
    img = self.readData(fn, z = (3,8));  
    print "Image size: " + str(img.shape)
    
    diff = img - data[:,:,3:8];
    print (diff.max(), diff.min())

    
    fn = os.path.join(basedir,'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z1000.ome.tif')        
    
    ds = self.dataSize(fn);
    print "Image size form dataSiZe: " + str(ds)    
    
    
    img = self.readData(fn);  
    print "Image size: " + str(img.shape)
    
    #dataSize
    print "dataSize  is %s" % str(self.dataSize(fn))
    print "dataZSize is %s" % str(self.dataZSize(fn))
    
    print "dataSize  is %s" % str(self.dataSize(fn, y = (10,20)))
    print "dataZSize is %s" % str(self.dataZSize(fn))
        
    #test writing multi channel image
    x = numpy.random.rand(50,100,30,4) * 10;
    x = x.astype('float32');
    self.writeData(os.path.join(basedir,'Test/Data/Tif/composite.tif'), x)

    y = self.readData(os.path.join(basedir,'Test/Data/Tif/composite.tif'));
    y.shape
        
    diff = x - y;
    print (diff.max(), diff.min())


    ### Meta data
    fn = os.path.join(basedir, 'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z1000.ome.tif')
    print self.readMetaData(fn, info = ['overlap', 'resolution'])

    fn = os.path.join(basedir, 'Test/Data/Tif/test.tif')
    print self.readMetaData(fn, info = ['overlap', 'resolution'])

    print self.readMetaData(fn)


if __name__ == "__main__":
    
    test();
    



# Test various tif libs:
#import numpy;
#import iDISCO.IO.IO as io
#x = numpy.random.rand(30,50,100,4) * 10;
#x = x.astype('float32');
#
#io.writeData('test.tif', x.astype('float32'))
#
#import tifffile as tiff
#
#tiff.imsave('test.tif', d)   
#
#t = tiff.TiffWriter('test.tif', bigtiff = True);
#
#t.save(x.astype('int16'), photometric = 'minisblack',  planarconfig = 'contig')
##t.save(x.astype('int16'), photometric = 'minisblack',  planarconfig = None)
#t.close();
#
#d = tiff.imread('Composite.tif')
#
#from libtiff import TIFF
#tif = TIFF.open('Composite.tif', mode='r')
#image = tif.read_image()
#image.shape
#
#from PIL import Image
#im = Image.open('Composite.tif');
#img = numpy.array(im);
#img.shape
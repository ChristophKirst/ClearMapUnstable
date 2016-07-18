# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 17:21:40 2016

@author: ckirst
"""

import os

import ClearMap.Settings as settings
import ClearMap.IO as io


fn = os.path.join(settings.ClearMapPath, 'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z1000.ome.tif')

io.readMetaData(fn, info = ['overlap', 'resolution', 'size'])

fn = os.path.join(settings.ClearMapPath, 'Test/Data/Tif/test.tif')

io.readMetaData(fn, info = ['overlap', 'resolution', 'size', 'description'])

io.readMetaData(fn)


fn = os.path.join(settings.ClearMapPath, 'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z\d{4}.ome.tif')

io.readMetaData(fn, info = ['overlap', 'resolution', 'size'])


import tifffile as t;


fn = os.path.join(settings.ClearMapPath, 'Test/Data/OME/16-17-27_0_8X-s3-20HF_UltraII_C00_xyz-Table Z1000.ome.tif')

tf = t.TiffFile(fn)

p = tf.pages[0];
o = p.tags['image_description']
xmlstring = o.value;

from lxml import etree

xml = etree.fromstring(xmlstring);

e1 = [x for x in xml.iter('{*}xyz-Table_X_Overlap')];
e2 = [x for x in xml.iter('{*}xyz-Table_Y_Overlap')];

e10 = e1[0];
e10.attrib['Value'] = '10.0';

xmlstring2 = etree.tostring(xml, pretty_print = False);



fnout = 'test.tif';

import ClearMap.IO as io
data = io.readData(fn);
tout = t.TiffWriter(fnout);
tout.save(data, description=xmlstring2);



tn = t.TiffFile(fnout)

p = tn.pages[0];
o = p.tags['image_description']
xmlstring = o.value;

xml = etree.fromstring(xmlstring);

e0 = [x for x in xml.iter('{*}xyz-Table_X_Overlap')];
e0[0].attrib


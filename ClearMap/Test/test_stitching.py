# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 17:11:51 2016

@author: ckirst
"""

import numpy as np

import ClearMap.Alignment.Stitching as st


imagepath = '/home/ckirst/Science/Projects/BrainActivityMap/Experiment/vasculature/160325_mosaic_13-43-57/';

imagename = os.path.join(imagepath, '13-43-57_mosaic_UltraII[00 x 00]_C00_xyz-Table Z0000.ome.tif')


d = st.findFileInfo(imagename);
print d


imageexpr = os.path.join(imagepath, '13-43-57_mosaic_UltraII\[(?P<col>\d{2}) x (?P<row>\d{2})\]_C00_xyz-Table Z(?P<z>\d{4}).ome.tif')

fns, infos = st.findFileList(imageexpr, groups = ('col', 'row', 'z'))
infos[0][1]= None

infos = np.array(infos)


b = np.zeros(infos.shape[0], dtype = bool);
for i in range(infos.shape[1]):
  b = np.logical_or(b, np.equal(infos[:,i], None))
idx = np.nonzero(b)


  

# test 
idx = np.logical_or(np.logical_or(np.equal(infos[:,0], None), np.equal(infos[:,1], None)), 
 np.nonzero( np.equal(infos[:,1], None))
















import os
import xml.etree.ElementTree as ET

from ClearMap.Settings import ClearMapPath


fn = os.path.join(ClearMapPath, 'Data/xml_import.xml');


tree = ET.parse(fn)
root = tree.getroot()


# read image info / tile info

imagepath = '/home/ckirst/Science/Projects/BrainActivityMap/Experiment/vasculature/160325_mosaic_13-43-57/';

imagepath = '/run/media/ckirst/New Volume/superTH/195A/160226_TH-s3-ov45-na004_14-22-21/'


imagename = os.path.join(imagepath, '13-43-57_mosaic_UltraII[00 x 00]_C00_xyz-Table Z0000.ome.tif')

imagename = os.path.join(imagepath, '14-22-21_TH-s3-ov45-na004_UltraII[00 x 00]_C00_xyz-Table Z0561.ome.tif')

import tifffile as tf

t = tf.TiffFile(imagename)

terastitcher --import --volin=<my_data> --volin_plugin="TiledXY|2Dseries" --imin_plugin="tiff2D" --ref1=y --ref2=-x --ref3=z --vxl1=0.8 --vxl2=0.8 --vxl3=1 --projout=xml_import


from PIL import Image
im = Image.open(imagename)


import re

fn = '13-43-57_mosaic_UltraII[00 x 00]_C00_xyz-Table Z0000.ome.tif'





s = re.compile(r'.*\[(?P<row>\d{2}) x (?P<col>\d{2})\].*(?P<z>\d{4}).*').search
m = s(fn); 

m.groupdict()




st.firstFile(fe)

fn = '13-43-57_mosaic_UltraII[00 x 00]_C00_xyz-Table Z0000.ome.tif'
fe = '13-43-57_mosaic_UltraII\[\d{2} x \d{2}\]_C00_xyz-Table Z\d{4}.ome.tif'

s = re.compile(fe).search
m = s(fn); 








import PIL.Image
from PIL.ExifTags import TAGS
import os

imagepath = '/home/ckirst/Science/Projects/BrainActivityMap/Experiment/vasculature/160325_mosaic_13-43-57/';
imagename = os.path.join(imagepath, '13-43-57_mosaic_UltraII[00 x 00]_C00_xyz-Table Z0000.ome.tif')

img = PIL.Image.open(imagename)

d = img.tag.as_dict();

for k in d.keys():
  print TAGS.get(k)
  
  

import xml.dom.minidom

xmls = xml.dom.minidom.parseString(d[270][0]) # or xml.dom.minidom.parseString(xml_string)
pretty_xml_as_string = xmls.toprettyxml()



f = open('test.txt', 'w');
f.write(pretty_xml_as_string);
f.close();




import re

xmlstring = d[270][0]; # or xml.dom.minidom.parseString(xml_string)
#xmlstring = re.sub(' xmlns="[^"]+"', '', xmlstring, count=1)
#xmlstring = re.sub(' xmlns:ca="[^"]+"', '', xmlstring, count=1)
#xmlstring = re.sub(' xmlns:xsi=="[^"]+"', '', xmlstring, count=1)


import xml.etree.ElementTree as et

root = et.fromstring(xmlstring);
 
 
 

from cStringIO import StringIO
from xml.etree import ElementTree

xmlin = """
<document xmlns="http://www.lol.com/"
          xmlns:foo="http://www.foo.com/"
          xmlns:bar="http://www.bar.com/">
    <something>
        <foo:thing>Hello,</foo:thing>
        <bar:thing>world!</bar:thing>
    </something>
</document>
"""

xml = None
namespaces = {}

for event, elem in ElementTree.iterparse(StringIO(xmlstring), ('start', 'start-ns')):
    if event == 'start-ns':
        if elem[0] in namespaces and namespaces[elem[0]] != elem[1]:
            # NOTE: It is perfectly valid to have the same prefix refer
            #     to different URI namespaces in different parts of the
            #     document. This exception serves as a reminder that this
            #     solution is not robust.    Use at your own peril.
            raise KeyError("Duplicate prefix with different URI found.")

        namespaces[str(elem[0])] = elem[1]

    elif event == 'start':
        if xml is None:
            xml = elem
            break

for prefix, uri in namespaces.iteritems():
    ElementTree.register_namespace(prefix, uri)

print ElementTree.tostring(xml)

namespaces = {None: 'http://www.openmicroscopy.org/Schemas/OME/2008-02', 'ca': 'http://www.openmicroscopy.org/Schemas/CA/2008-02', 'xsi': 'http://www.w3.org/2001/XMLSchema-instance'}

root.findall('xyz-Table_X_Overlap', namespaces = namespaces)



from lxml import etree

xmlt = etree.fromstring(str(xmlstring))

for el in xmlt.iter('{*}xyz-Table_X_Overlap'):


e1 = [x for x in xmlt.iter('{*}Pixels')];
e1[0].attrib





import os

imagepath = '/home/ckirst/Science/Projects/BrainActivityMap/Experiment/vasculature/160325_mosaic_13-43-57/';
imagename = os.path.join(imagepath, '13-43-57_mosaic_UltraII[00 x 00]_C00_xyz-Table Z0000.ome.tif')


import ClearMap.Alignment.Stitching as st;
reload(st);

d = st.findFileInfo(imagename);
print d


xmlstring = """<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE TeraStitcher SYSTEM "TeraStitcher.DTD">
<TeraStitcher volume_format="TiledXY|2Dseries">
    <stacks_dir value="C:\Campus BioMedico\Stitching\Data\TeraStitcher_TESTSET\tomo300511_subv2.sparse3D" />
    <voxel_dims V="0.8" H="-0.8" D="1" />
    <origin V="7.2" H="15.8088" D="13.899" />
    <mechanical_displacements V="300" H="-300" />
    <dimensions stack_rows="2" stack_columns="2" stack_slices="600" />
    <STACKS>
        <Stack ROW="0" COL="0" ABS_V="0" ABS_H="0" ABS_D="0" STITCHABLE="no" DIR_NAME="072000/072000_154000" IMG_REGEX="[0-9]+_[0-9]+_[0-9]+.tiff" Z_RANGES="[0,600)">
            <NORTH_displacements />
            <EAST_displacements />
            <SOUTH_displacements />
            <WEST_displacements />
        </Stack>
        <Stack ROW="0" COL="1" ABS_V="0" ABS_H="375" ABS_D="0" STITCHABLE="no" DIR_NAME="072000/072000_151000" IMG_REGEX="[0-9]+_[0-9]+_[0-9]+.tiff" Z_RANGES="[0,600)">
            <NORTH_displacements />
            <EAST_displacements />
            <SOUTH_displacements />
            <WEST_displacements />
        </Stack>
        <Stack ROW="1" COL="0" ABS_V="375" ABS_H="0" ABS_D="0" STITCHABLE="no" DIR_NAME="075000/075000_154000" IMG_REGEX="[0-9]+_[0-9]+_[0-9]+.tiff" Z_RANGES="[9,135);[162,207);[225,585)">
            <NORTH_displacements />
            <EAST_displacements />
            <SOUTH_displacements />
            <WEST_displacements />
        </Stack>
        <Stack ROW="1" COL="1" ABS_V="375" ABS_H="375" ABS_D="0" STITCHABLE="no" DIR_NAME="075000/075000_151000" IMG_REGEX="[0-9]+_[0-9]+_[0-9]+.tiff" Z_RANGES="[0,600)">
            <NORTH_displacements />
            <EAST_displacements />
            <SOUTH_displacements />
            <WEST_displacements />
        </Stack>
    </STACKS>
</TeraStitcher>
"""

from lxml import etree

xmlt = etree.fromstring(str(xmlstring))

print(etree.tostring(xmlt, pretty_print=False, xml_declaration=True, encoding="utf-8",  doctype='<!DOCTYPE TeraStitcher SYSTEM "TeraStitcher.DTD">'))


# create a stack entry







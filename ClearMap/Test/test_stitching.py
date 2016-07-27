# -*- coding: utf-8 -*-
"""
Created on Fri Jul  8 17:11:51 2016
"""    
    
    
import ClearMap.Alignment.Stitching as st
reload(st)

st.printSettings()


import os


datadir = '/home/mtllab/Documents/th/cropped';
#datadir = '/home/mtllab/Documents/th/160412_mosaic_15-20-19';
outdir = '/home/mtllab/Documents/th/stitched';

#import glob
#fn = os.path.join(datadir, '14-22-21_TH-s3-ov45-na004_UltraII[[]* x *[]]_C00_xyz-Table Z*.ome.tif')
#fn = os.path.join(datadir, '*')
#print fn
#print glob.glob(fn)

#fileExpression = os.path.join(datadir, r'15-20-19_mosaic_UltraII\[(?P<row>\d{2}) x (?P<col>\d{2})\]_C00_xyz-Table Z(?P<z>\d{4}).ome.tif')
fileExpression = os.path.join(datadir, r'15-20-19_mosaic_UltraII_(?P<row>\d{2})_x_(?P<col>\d{2})_C00_xyz-Table Z(?P<z>\d{4}).ome.tif')

#st.findFileList(fileExpression)

ff = st.findFirstFile(fileExpression, sort = True)

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
alignfile = st.alignData(importfile, search = (75,75,0), xmlResultFile = os.path.join(datadir, 'TeraStitcher_align.xml'), subRegion=((all, all),(all,all),(140, 150)));

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
result = st.stitchData(placefile, resultPath = os.path.join(datadir, 'stiched_noblend.tif'), bitDepth=16, algorithm='NOBLEND')
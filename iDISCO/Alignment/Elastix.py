# -*- coding: utf-8 -*-

"""
Interface to Elastix for brain alignment

Notes:
    - elastix finds map T: fixed -> moving
    - the result folder may contain an image (mhd file) that is T^-1(moving) i.e. has size of fixed image
    - transformix applied to data gives T^(-1)(data) !
    - transformix applied to points gives T(points)
    - point arrays are assumed to be in (x,y,z) coordinates consistent with (x,y,z) array represenation of images in iDISCO

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import os
import tempfile
import shutil
import numpy
import re

import iDISCO.Settings as settings

import iDISCO.IO.IO as io

##############################################################################
# Initialization and Enviroment Settings
##############################################################################

ElastixBinary = None;
ElastixLib = None;
TransformixBinary = None;
    
Initialized = False;
    
def printSettings():
    global ElastixBinary, ElastixLib, TransformixBinary, Initialized
    
    if Initialized:
        print "ElastixBinary     = %s" % ElastixBinary;
        print "ElastixLib        = %s" % ElastixLib;
        print "TransformixBinary = %s" % TransformixBinary;
    else:
        print "Elastix not initialized";


def setElastixLibraryPath(path = None): 
    """Add elastix libs to path in linux"""
     
    if path is None:
        path = settings.ElastixPath;
    
    if os.environ.has_key('LD_LIBRARY_PATH'):
        lp = os.environ['LD_LIBRARY_PATH'];
        if not path in lp.split(':'):
            os.environ['LD_LIBRARY_PATH'] = lp + ':' + path;
    else:
        os.environ['LD_LIBRARY_PATH'] = path


def initializeElastix(path = None):
    """Initialize all paths and binaries of elastix using installation directory elastixdir"""
    
    global ElastixBinary, ElastixLib, TransformixBinary, Initialized
    
    if path is None:
        path = settings.ElastixPath;
    
    #search for elastix binary
    elastixbin = os.path.join(path, 'bin/elastix');
    if os.path.exists(elastixbin):
        ElastixBinary = elastixbin;
    else:
        raise RuntimeError("Cannot find elastix binary %s, set path in Settings.py accordingly!" % elastixbin);
    
    #search for transformix binarx
    transformixbin = os.path.join(path, 'bin/transformix');
    if os.path.exists(transformixbin):
        TransformixBinary = transformixbin;
    else:
        raise RuntimeError("Cannot find transformix binary %s set path in Settings.py accordingly!" % transformixbin);
    
    #search for elastix libs
    elastixlib = os.path.join(path, 'lib');
    if os.path.exists(elastixlib):
        ElastixLib = elastixlib;
    else:
        raise RuntimeError("Cannot find elastix libs in %s  set path in Settings.py accordingly!" % elastixlib);
    
    #set path
    self.setElastixLibraryPath(elastixlib);
        
    Initialized = True;
    
    print "Elastix sucessfully initialized from path: %s" % path;
    
    return path;



initializeElastix();


def checkElastixInitialized():
    """Checks if elastix is initialized"""
    
    global Initialized;
    
    if not Initialized:
        raise RuntimeError("Elastix not initialized: run initializeElastix(path) with proper path to elastix first");
    #print ElastixSettings.ElastixBinary;

    return True;



##############################################################################
# Basic interface routines
##############################################################################

def getTransformParameterFile(resultdir):
    """Finds and returns the transformation parameter file generated by elastix"""    
    
    files = os.listdir(resultdir);
    files = [x for x in files if re.match('TransformParameters.\d.txt', x)];
    files.sort();
    
    if files == []:
        raise RuntimeError('Cannot find a valid transformation file in ' + resultdir);
    
    return os.path.join(resultdir, files[-1])


def setPathTransformParameterFiles(resultdir):
    """Replaces relative with abolsute path in the parameter files in the result directory"""
    print resultdir 
    files = os.listdir(resultdir);
    files = [x for x in files if re.match('TransformParameters.\d.txt', x)];
    files.sort();
    
    if files == []:
        raise RuntimeError('Cannot find a valid transformation file in ' + resultdir);
    
    rec = re.compile("\(InitialTransformParametersFileName \"(?P<parname>.*)\"\)");
    
    for f in files:
        fh, tmpfn = tempfile.mkstemp();
        ff = os.path.join(resultdir, f);
        #print ff        
        
        with open(tmpfn, 'w') as newfile:
            with open(ff) as parfile:
                for line in parfile:
                    #print line
                    m = rec.match(line);
                    if m != None:
                        pn = m.group('parname');
                        if pn != 'NoInitialTransform':
                            pathn, filen = os.path.split(pn);
                            filen = os.path.join(resultdir, filen);
                            newfile.write(line.replace(pn, filen));
                        else:
                            newfile.write(line);
                    else:
                        newfile.write(line);
                            
        os.close(fh);
        os.remove(ff);
        shutil.move(tmpfn, ff);


def parseElastixOutputPoints(filename, indices = True):
    """Parses the output points from the output file of transformix"""
    
    with open(filename) as f:
        lines = f.readlines()
        f.close();
    
    np = len(lines);
    
    if np == 0:
        return numpy.zeros((0,3));
    
    points = numpy.zeros((np, 3));
    k = 0;
    for line in lines:
        ls = line.split();
        if indices:
            for i in range(0,3):
                points[k,i] = float(ls[i+22]);
        else:
            for i in range(0,3):
                points[k,i] = float(ls[i+30]);
        
        k += 1;
    
    return points;
          
         
def getTransformFileSizeAndSpacing(transformfile):
    """Parse the image size and spacing from a transformation parameter file"""
    
    resi = re.compile("\(Size (?P<size>.*)\)");
    resp = re.compile("\(Spacing (?P<spacing>.*)\)");
    
    si = None;
    sp = None;
    with open(transformfile) as parfile:
        for line in parfile:
            #print line;
            m = resi.match(line)
            if m != None:
                pn = m.group('size');
                si = pn.split();
                #print si
                
            m = resp.match(line);
            if m != None:
                pn = m.group('spacing');
                sp = pn.split();
                #print sp 
    
        parfile.close();
    
    si = [float(x) for x in si];
    sp = [float(x) for x in sp];
    
    return si, sp


def getResultDataFile(resultdir):
    """Returns the mhd result file in a result directory"""
    
    files = os.listdir(resultdir);
    files = [x for x in files if re.match('.*.mhd', x)];
    files.sort();
    
    if files == []:
        raise RuntimeError('Cannot find a valid result data file in ' + resultdir);
    
    return os.path.join(resultdir, files[0])


    
def setTransformFileSizeAndSpacing(transformfile, size, spacing):
    """Replaces size and scale in the transformation parameter file"""
    
    resi = re.compile("\(Size (?P<size>.*)\)");
    resp = re.compile("\(Spacing (?P<spacing>.*)\)");
    
    fh, tmpfn = tempfile.mkstemp();
    
    si = [int(x) for x in size];
    
    with open(transformfile) as parfile:        
        with open(tmpfn, 'w') as newfile:
            for line in parfile:
                #print line
                
                m = resi.match(line)
                if m != None:
                    newfile.write("(Size %d %d %d)" % si);
                else:
                    m = resp.match(line)
                    if m != None:
                        newfile.write("(Spacing %d %d %d)" % spacing);
                    else:
                        newfile.write(line);
            
            newfile.close();               
            parfile.close();
            
            os.remove(transformfile);
            shutil.move(tmpfn, transformfile);
        


def rescaleSizeAndSpacing(size, spacing, scale):
    si = [int(x * scale) for x in size];
    sp = spacing / scale;
    
    return si, sp



##############################################################################
# Elastix Runs
##############################################################################

def alignData(fixedImage, movingImage, affineParameterFile, bSplineParameterFile = None, resultDirectory = None):
    """Align images using elastix, estimates transformation T:fixed -> moving"""
    
    self.checkElastixInitialized();
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
    """Transform a raw data set to reference using the elastix alignment results"""
    
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
        transformParameterFile = self.getTransformParameterFile(transformparameterdir);
    else:
        transformparameterdir = os.path.split(transformParameterFile);
        transformparameterdir = transformparameterdir[0];
    
    #transform
    #make path in parameterfiles absolute
    self.setPathTransformParameterFiles(transformparameterdir);
   
    #transformix -in inputImage.ext -out outputDirectory -tp TransformParameters.txt
    cmd = TransformixBinary + ' -in ' + imgname + ' -out ' + resultdirname + ' -tp ' + transformParameterFile;
    
    res = os.system(cmd);
    
    if res != 0:
        raise RuntimeError('transformData: failed executing: ' + cmd);
    
    
    if not isinstance(source, basestring):
        os.remove(imgname);

    if sink == []:
        return self.getResultDataFile(resultdirname);
    elif sink is None:
        resultfile = self.getResultDataFile(resultdirname);
        return io.readData(resultfile);
    elif isinstance(sink, basestring):
        resultfile = self.getResultDataFile(resultdirname);
        return io.convertData(resultfile, sink);
    else:
        raise RuntimeError('transformData: sink not valid!');


def deformationField(sink = [], transformParameterFile = None, transformDirectory = None, resultDirectory = None):
    """Create the deformation field T(x) - x"""
    
    global TransformixBinary;    
    
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
        transformParameterFile = self.getTransformParameterFile(transformparameterdir);
    else:
        transformparameterdir = os.path.split(transformParameterFile);
        transformparameterdir = transformparameterdir[0];
    
    #transform
    #make path in parameterfiles absolute
    self.setPathTransformParameterFiles(transformparameterdir);
   
    #transformix -in inputImage.ext -out outputDirectory -tp TransformParameters.txt
    cmd = TransformixBinary + ' -def all -out ' + resultdirname + ' -tp ' + transformParameterFile;
    
    res = os.system(cmd);
    
    if res != 0:
        raise RuntimeError('deformationField: failed executing: ' + cmd);
    
    
    if sink == []:
        return self.getResultDataFile(resultdirname);
    elif sink is None:
        resultfile = self.getResultDataFile(resultdirname);
        data = io.readData(resultfile);
        if resultDirectory is None:
            shutil.rmtree(resultdirname);
        return data;
    elif isinstance(sink, basestring):
        resultfile = self.getResultDataFile(resultdirname);
        data = io.convertData(resultfile, sink);
        if resultDirectory is None:
            shutil.rmtree(resultdirname);
        return data;
    else:
        raise RuntimeError('deformationField: sink not valid!');


def deformationDistance(deformationField, sink = None, scale = None):
    """Compute the distance field from a deformation vector field"""
    
    df = numpy.square(deformationField);
    if not scale is None:
        for i in range(3):
            df[:,:,:,i] = df[:,:,:,i] * (scale[i] * scale[i]);
    return numpy.sqrt(numpy.sum(df, axis = 3));
    

def writePoints(filename, points, indices = True):
    """Write points as elastix/transformix point file"""

    #points = points[:,[1,0,2]]; # points in iDISCO (y,x,z) -> permute to (x,y,z)
  
    with open(filename, 'w') as pointfile:
        if indices:
            pointfile.write('index\n')
        else:
            pointfile.write('point\n')
    
        pointfile.write(str(points.shape[0]) + '\n');
        numpy.savetxt(pointfile, points, delimiter = ' ', newline = '\n', fmt = '%.5e')
        pointfile.close();
    
    return filename;



def transformPoints(source, sink = None, transformParameterFile = None, transformDirectory = None, indices = True, resultDirectory = None, tmpFile = None):
    """Transform pixel coordinates of cell centers using a calculated transformation obtained form elastix"""
    """Note: assumes iDISCOs (y,x,z) representation of points as arrays/sources but (x,y,z) representation of elastix type point files""" 
    
    global TransformixBinary;    
    
    self.checkElastixInitialized();    
    global ElastixSettings;

    if tmpFile == None:
        tmpFile = os.path.join(tempfile.tempdir, 'elastix_input.txt');

    # write text file
    if isinstance(source, basestring):
        
        #check if we have elastix signature                 
        with open(source) as f:
            line = f.readline();
            f.close();
            
            if line[:5] == 'point' or line[:5] != 'index':
                txtfile = source;
            else:                
                points = io.readPoints(source);
                #points = points[:,[1,0,2]];
                txtfile = tmpFile;
                self.writePoints(txtfile, points); 
    
    elif isinstance(source, numpy.ndarray):
        txtfile = tmpFile;
        #points = source[:,[1,0,2]];
        self.writePoints(txtfile, source);
        
    else:
        raise RuntimeError('transformPoints: source not string or array!');
    
    
    if resultDirectory == None:
        outdirname = os.path.join(tempfile.tempdir, 'elastix_output');
    else:
        outdirname = resultDirectory;
        
    if not os.path.exists(outdirname):
        os.makedirs(outdirname);
        
    
    if transformParameterFile == None:
        if transformDirectory == None:
            RuntimeError('neither alignment directory and transformation parameter file specified!'); 
        transformparameterdir = transformDirectory
        transformparameterfile = getTransformParameterFile(transformparameterdir);
    else:
        transformparameterdir = os.path.split(transformParameterFile);
        transformparameterdir  = transformparameterdir[0];
        transformparameterfile = transformParameterFile;
    
    #transform
    #make path in parameterfiles absolute
    self.setPathTransformParameterFiles(transformparameterdir);
    
    #run transformix   
    cmd = TransformixBinary + ' -def ' + txtfile + ' -out ' + outdirname + ' -tp ' + transformparameterfile;
    res = os.system(cmd);
    
    if res != 0:
        raise RuntimeError('failed executing ' + cmd);
    
    
    #read data / file 
    if sink == []:
        return io.path.join(outdirname, 'outputpoints.txt')
    
    else:
        #read coordinates
        transpoints = parseElastixOutputPoints(os.path.join(outdirname, 'outputpoints.txt'), indices = indices);

        #correct x,y,z to y,x,z
        #transpoints = transpoints[:,[1,0,2]];     
        
        #cleanup
        for f in os.listdir(outdirname):
            os.remove(os.path.join(outdirname, f));
        os.rmdir(outdirname)
        
        return io.writePoints(sink, transpoints);

        
        
 

##############################################################################
# Test
##############################################################################

     
def test():
    """Test Elastix module"""
    import iDISCO.Alignment.Elastix as self
    reload(self)
    
    from iDISCO.Settings import IDISCOPath;
    import os, numpy
    
    p = IDISCOPath;
    
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
    df =  '/tmp/elastix_output/deformationField.mhd';

    import iDISCO.IO.IO as io
    data = io.readData('/tmp/elastix_output/deformationField.mhd');
    
    ds = self.deformationDistance(data);
    
    io.writeData(os.path.join(p, 'Test/Elastix/Output/distances.raw'), ds);
    
if __name__ == "__main__":
    self.test();
    





  

    
    
    

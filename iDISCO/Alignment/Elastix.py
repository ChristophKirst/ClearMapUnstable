# -*- coding: utf-8 -*-

"""
Interface to Elastix for brain alignment

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

from iDISCO.Parameter import AlignmentParameter


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
    
    files = os.listdir(resultdir);
    files = [x for x in files if re.match('TransformParameters.\d.txt', x)];
    files.sort();
    
    if files == []:
        raise RuntimeError('Cannot find a valid transformation file in ' + resultdir);
    
    rec = re.compile("\(InitialTransformParametersFileName \"(?P<parname>.*)\"\)");
    
    for f in files:
        fh, tmpfn = tempfile.mkstemp();
        ff = os.path.join(resultdir, f);
        
        print ff        
        
        with open(tmpfn, 'w') as newfile:
            with open(ff) as parfile:
                for line in parfile:
                    print line
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
    
    print np
    
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
          
          
def setElastixLibraryPath(parameter = AlignmentParameter()): 
    #make libs for elastix visible in linux
    lp = os.environ['LD_LIBRARY_PATH']   
    os.environ['LD_LIBRARY_PATH'] = lp + ':' + parameter.Lib;
    


def transformPoints(points, transformparameterfile = None, parameter = AlignmentParameter(), read = True, tmpfile = None, outdir = None, indices = True):
    """Transform pixel coordinates of cell centers using a calculated transformation obtained fomr elastix"""
    
    # write text file
    if isinstance(points, basestring):
        txtfile = points;
        #check if we have elastix signature
        
        with open(txtfile) as f:
            line = f.readline();
            f.close();
            
            if line[:5] != 'point' and line[:5] != 'index':
                fh, tmpfn = tempfile.mkstemp();
                with open(txtfile) as pointfile:
                    with open(tmpfn, 'w') as newfile:
                        lines = pointfile.readlines();
                        
                        if indices:
                            newfile.write('index\n')
                        else:
                            newfile.write('point\n')
                        
                        newfile.write(str(lines.shape[0]) + '\n');
                        newfile.writelines(lines);
                
                txtfile = tmpfn;   
        
    else: # array of points
        if tmpfile == None:
            txtfile = os.path.join(tempfile.tempdir, 'elastix_input.txt');
        else:
            txtfile = tempfile;
        
        with open(txtfile, 'w') as pointfile:
            if indices:
                pointfile.write('index\n')
            else:
                pointfile.write('point\n')
        
            pointfile.write(str(points.shape[0]) + '\n');
            numpy.savetxt(pointfile, points, delimiter = ' ', newline = '\n', fmt = '%.5e')
            pointfile.close();

    
    if outdir == None:
        outdirname = os.path.join(tempfile.tempdir, 'elastix_output');
    else:
        outdirname = outdir;
        
    if not os.path.exists(outdirname):
        os.makedirs(outdirname);
        
    
    if transformparameterfile == None:
        transformparameterdir = parameter.AlignmentDirectory
        transformparameterfile = getTransformParameterFile(transformparameterdir);
    else:
        transformparameterdir = os.path.split(transformparameterfile);
    
    #transform
    #make path in parameterfiles absolute
    self.setPathTransformParameterFiles(transformparameterdir);
    
    #set library path
    self.setElastixLibraryPath(parameter);

    #run transformix   
    cmd = parameter.Transformix + ' -def ' + txtfile + ' -out ' + outdirname + ' -tp ' + transformparameterfile;
    res = os.system(cmd);
    
    if res != 0:
        raise RuntimeError('failed executing ' + cmd);
    
    
    #read data / file 
    if outdir == None:    
        #read coordinates
        transpoints = parseElastixOutputPoints(os.path.join(outdirname, 'outputpoints.txt'), indices = indices);
        
        #cleanup
        for f in os.listdir(outdirname):
            os.remove(os.path.join(outdirname, f));
        os.rmdir(outdirname)
        
        return transpoints;
    
    else:
        return outdirname + 'outputpoints.txt'


def downSampleData(filename):
    """Samples the raw data to match the reference data for alignmnet"""
    
    print 'todo';
 

def alignData(movingimage, fixedimage, parameter = AlignmentParameter()):
    """Align data to reference using elastix"""
    
    affineparfile = parameter.AffineParameterFile;
    bsplinefile   = parameter.BSplineParameterFile;

    outdir = parameter.AlignmentDirectory;
    
    if outdir == None:
        outdir = tempfile.gettempdir();
    
    if not os.path.exists(outdir):
        os.mkdir(outdir);
    
    cmd = parameter.Elastix + ' -threads 16 -m ' + movingimage + ' -f ' + fixedimage + ' -p ' + affineparfile + ' -p ' + bsplinefile + ' -out ' + outdir;
    #$ELASTIX -threads 16 -m $MOVINGIMAGE -f $FIXEDIMAGE -fMask $FIXEDIMAGE_MASK -p  $AFFINEPARFILE -p $BSPLINEPARFILE -out $ELASTIX_OUTPUT_DIR
    
    res = os.system(cmd);
    
    if res != 0:
        raise RuntimeError('failed executing ' + cmd);
    
    return outdir


     
def test():
    from iDISCO.Parameter import iDISCOPath;
    p = iDISCOPath();
    
    resultdir = os.path.join(p, 'Test/Elastix/Output');
    
    print 'Searching for transformation parameter file in ' + resultdir;
    pf = self.getTransformParameterFile(resultdir)
      
    print 'Found: ' + pf;
    
    
    #replace path in trasform parameter files:
    self.setPathTransformParameterFiles(resultdir)
    
    #transform
    pts = numpy.random.rand(5,3);    
    parameter = AlignmentParameter();
    parameter.Lib = '/home/ckirst/programs/elastix/lib';
    parameter.Transformix = '/home/ckirst/programs/elastix/bin/transformix';
    
    print 'Transforming points: '
    tpts = self.transfromByEalstix(pts, pf, indices = False, parameter= parameter);
    print pts
    print 'Transformed points: '
    print tpts
    
    
if __name__ == "__main__":
    self.test();
    





  

    
    
    

# -*- coding: utf-8 -*-
"""
Simple Interface to Imaris Files 

read data and write points to Imarise files

Note: to write points make sure the original file has at least one spot object

Todo: fix writing new spots to imaris file
      get settings to directly render points as 'pixel' and not as spheres

Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import h5py
import numpy

import sys
self = sys.modules[__name__];

import iDISCO.IO.IO as io


def openFile(filename, mode = "a"):
    """Open Imaris file as hdf5 object"""
    
    return h5py.File(filename, mode);

    
def closeFile(h5file, mode = "a"):
    """Close Imaris hdf5 file object"""
    
    return h5file.close();


def readDataSet(h5file, resolution = 0, channel = 0, timepoint = 0):
    """Open Imaris file and returns hdf5 dataset object"""
    
    dsname = "/DataSet/ResolutionLevel " + str(resolution) + "/TimePoint " + str(timepoint) + "/Channel " + str(channel) + "/Data";
    return h5file.get(dsname);


def dataSize(filename, resolution = 0, channel = 0, timepoint = 0, **args):
    f = self.openFile(filename);
    
    #todo: read array size not full data !
    ds = self.readDataSet(f, resolution = 0, channel = 0, timepoint = 0);
    dims = list(ds.shape);
    #dims = (dims[1], dims[2], dims[0]);
    
    return io.dataSizeFromDataRange(dims, **args);

 
def dataZSize(filename, **args):
    dsize = self.dataSize(filename, **args);
    return dsize[2];


def readData(filename, x = all, y = all, z = all, resolution = 0, channel = 0, timepoint = 0, **args):
    
    f = h5py.File(filename, "r");
    dataset = self.readDataSet(f, resolution = resolution, channel = channel, timepoint = timepoint);
    dsize = dataset.shape;
    
    rz = io.toDataRange(dsize[0], r = z);
    ry = io.toDataRange(dsize[1], r = y);
    rx = io.toDataRange(dsize[2], r = x);    
    
    data = dataset[rz[0]:rz[1],ry[0]:ry[1],rx[0]:rx[1]];
    data = data.transpose((2,1,0)); # imaris stores files in reverse x,y,z ordering
    #data = dataset[x[0]:x[1],y[0]:y[1],z[0]:z[1]];
    
    f.close();
    
    return data;

   

def transformToImaris(points, scale = (4.0625, 4.0625, 3)):
    """Transform pixel coordinates of cell centers to work in Imaris"""
    if len(scale) == 1:
        scale = (scale, scale, scale);    
    
    for i in range(3):
        points[:,i] = points[:,i] * scale[i];
        
    return points


def writePoints(filename, points, mode = "o", radius = 0.5):
    """Write points to Imaris file"""
    
    if isinstance(filename, basestring):
        h5file = self.openFile(filename);
    else:
        h5file = filename;
      
    #delete Scene8 info so do not need to write it
    s8 = "/Scene8";
    if h5file.get(s8) != None:
        del h5file[s8];
    
    #get number of point sets
    ct = "Scene/Content";
    ca = h5file.get(ct).attrs
    np = ca.get("NumberOfPoints");
    if np == None:
        ca.create("NumberOfPoints", 1, dtype = "uint64")
        np = 0;
        mode = "a";
    
    #write points
    if mode == "a":  # add points -> need to test further
        
        #update number of points
        np = np + 1;
        ca["NumberOfPoints"] = np;
        
        # Create Pointset and attributes
        pn = ct + "/Points" + str(int(np-1));
        
        histlog = numpy.array([''], dtype='|S1');        
        name = numpy.array(['Spots ' + str(int(np))], dtype='|S' + str(6 + len(str(int(np)))));
        material = numpy.array([ '<bpMaterial mDiffusion="0.8 0.8 0.8" mSpecular="0.2 0.2 0.2" mEmission="0 0 0" mAmbient="0 0 0" mTransparency="0" mShinyness="0.1"/>\n'], dtype='|S133');       
        iid = 200000 + np # int64
        unit = numpy.array(['um'], dtype='|S2');
        descrp = numpy.array([''], dtype='|S1');
                
        pg = h5file.create_group(pn);
        pa = pg.attrs;
        pa.create("HistoryLog", histlog);
        pa.create("Name", name);        
        pa.create("Material", material);            
        pa.create("Id", iid, dtype = 'uint64');
        pa.create("Unit", unit);
        pa.create("Description", descrp);     
        
        #add TimeInfos
        
        tb = h5file.get('DataSetTimes/TimeBegin');
        tb = tb[0][1];
        
        #del h5file[pn + "/TimeInfos"]
        h5file.create_dataset(pn + "/TimeInfos", data = numpy.array([tb[0:23]], dtype='|S23'))
        
    else:
        pn = ct + "/Points" + str(int(np-1));
        
    
    #  make points
    npts = points.shape[0];
    
    if points.shape[1] != 3:
        raise StandardError("Points shape is not (n,3)!");
    
    points = points[:,[1,0,2]]; # todo: check exchange of coordinates
    pts = numpy.c_[points, radius * numpy.ones(npts)];
    ts =  numpy.zeros(npts);
    
    # write points
    pnt = pn + '/Time';
    pnc = pn + '/CoordsXYZR';
    
    if h5file.get(pnt) != None:
        del h5file[pnt];
    h5file.create_dataset(pnt, shape=(npts,1), dtype='i64', data=ts);
    
    
    if h5file.get(pnc) != None:
        del h5file[pnc];
    h5file.create_dataset(pnc, shape=pts.shape, dtype='f32', data=pts);
    
    if isinstance(filename, basestring):
        h5file.close();


def writeData(filename, **args):
    raise RuntimeError("Writing data to imaris file not implemented yet""")


def readPoints(filename):
    raise RuntimeError("Reading points form imaris file not implemented yet""")


def copyData(source, sink):
    io.copyFile(source, sink);


def test():
    """Test Imaris module"""
    import iDISCO.IO.IO as io
    import iDISCO.IO.Imaris as self
    
    from iDISCO.Parameter import iDISCOPath
    import os
    import numpy

    basedir = iDISCOPath();
    fn = os.path.join(basedir,'Test/Data/Imaris/test for spots with spot.ims')  
    fn = os.path.join(basedir,'Test/Data/Imaris/test for spots added spot.ims') 
    

    import h5py
    f = h5py.File(fn, "a");    
    
    dsname = "/DataSet/ResolutionLevel 0/TimePoint 0/Channel 0/Data"
    ds = f.get(dsname)


    data = numpy.random.rand(20,50,10);
    data[5:15, 20:45, 2:9] = 0;
    data = 20 * data;
    data = data.astype('int32');

    #reload(self)
    print "writing raw image to: " + fn;    
    self.writeData(fn, data);

    
    
    
if __name__ == "__main__":
    self.test();
    







"""


    import h5py    
    
    fn = "/home/ckirst/Data/Science/Projects/BrainActivityMap/iDISCO_2015_06/Adult cfos all shaved 20HF 150523.ims";
    
    dsname = "/DataSet/ResolutionLevel 4/TimePoint 0/Channel 0/Data"

    
    f = h5py.File(fn, "a");
    ds = f.get(dsname)


    fn = "/home/ckirst/Data/Science/Projects/BrainActivityMap/ImarisTest/test for spots added spot.ims";
    f = h5py.File(fn, "a");

### write spots to file with spots

#Scene
pn = '/Scene/Content/Points0'

pnt = pn + '/Time';
pnc = pn + '/CoordsXYZR';

pts = numpy.arange(1, 10);
npts = pts.size;

pts = numpy.array([1900 + 50 * numpy.sin(2 * numpy.pi * pts / npts), 2000 + 50 * numpy.cos(2 * numpy.pi * pts / npts), 1700  + numpy.zeros(npts), 2.5  + numpy.zeros(npts)])
pts = numpy.transpose(pts);
t = numpy.zeros(npts);

del f[pnt]
dt = f.create_dataset(pnt, shape=(npts,1), dtype='i64', data=t)

del f[pnc]
dc = f.create_dataset(pnc, shape=pts.shape, dtype='f32', data=pts)


fn = '/home/nicolas/Windows/Nico/cfosRegistrations/Adult cfos C row 20HF 150524.ims';

#
pts = np.arange(1, 10);
npts = pts.size;

pts2 = np.array([1900 + 50 * np.sin(2 * np.pi * pts / npts), 2000 + 50 * np.cos(2 * np.pi * pts / npts), 1700  + np.zeros(npts), 2.5  + np.zeros(npts)])
pts2 = np.transpose(pts2);
t = np.zeros(npts);

ti = '2015-03-13 21:16:16.000';
ti = np.array([ti.encode("ascii", "ignore")])

# open hdf5 file
f = h5py.File(fn, "a")
dt = f.create_dataset(pnt, shape=(npts,1), dtype='i64', data=t)
dc = f.create_dataset(pnc, shape=pts2.shape, dtype='f32', data=pts2)
di = f.create_dataset(pni, data = ti)

#attributes

f.close()

"""
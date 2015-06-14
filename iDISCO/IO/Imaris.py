# -*- coding: utf-8 -*-
"""
Simple Interface to Imaris Files 

read data and write points


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import h5py
import numpy

import sys
self = sys.modules[__name__];


def openFile(filename, mode = "a"):
    return h5py.File(filename, "a");
    
def closeFile(h5file, mode = "a"):
    return h5file.close();


def readDataSet(h5file, resolution = 4, channel = 0, timepoint = 0):
    """Read data from Imaris file and returns HDF5 dataset"""
    dsname = "/DataSet/ResolutionLevel " + str(resolution) + "/TimePoint " + str(timepoint) + "/Channel " + str(channel) + "/Data";
    return h5file.get(dsname);
    
    
def readData(filename, x = all, y = all, z = all, resolution = 4, channel = 0, timepoint = 0):
    
    f = h5py.File(filename, "r");
    dataset = self.readDataSet(f, resolution = resolution, channel = channel, timepoint  = timepoint);
    datasetsize = dataset.shape;
    
    if x == all:
        x = (0, datasetsize[0]);
    
    if y == all:
        y = (0, datasetsize[1]);
    
    if z == all:
        z = (0, datasetsize[2]);
    
    
    img = dataset[z[0]:z[1],y[0]:y[1],x[0]:x[1]];
    img = img.transpose((2,1,0)); # imaris stores files in reverse x,y,z ordering
    #img = dataset[1200:1400,1200:1400, zrange[0]:zrange[1]];   
    #img = dataset[0:50, 0:50, zrange[0]:zrange[1]];
    
    f.close();
    
    return img;



  
def writePoints(h5file, points, mode = "a", radius = 0.5):
    """Write points to Imaris file"""
    
    #delete Scene8 info to no need to write it
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

  
    
"""

if __name__ == "__main__":
    
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
    



# test

"""




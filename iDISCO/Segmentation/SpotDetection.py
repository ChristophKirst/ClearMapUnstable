# -*- coding: utf-8 -*-

"""
Spot Detection Algorithm

read data and write points


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import numpy
import math

#import mahotas as morph
import scipy.ndimage.morphology as morph
from scipy.ndimage.filters import correlate
#from scipy.ndimage import maximum_filter
from scipy.ndimage.measurements import label

from skimage.morphology import reconstruction
#from skimage.measure import regionprops
from scipy.ndimage.measurements import center_of_mass
from mahotas import regmax, regmin

import iDISCO.IO.Imaris as io


def structureElementOffsets(sesize):
    sesize = numpy.array(sesize);
    ndim = len(sesize);
    o = numpy.array(((0,0),(0,0), (0,0)));
    for i in range(0, ndim):   
        if sesize[i] % 2 == 0:
            o[i,0] = sesize[i]/2;
            o[i,1] = sesize[i]/2 - 1 + 1;
        else:
            o[i,0] = (sesize[i]-1)/2;
            o[i,1] = o[i,0] + 1;

    return o.astype('int');


def structureELement2D(setype = 'Disk', sesize = (3,3)):
    """Creates specific 2d structuring elements""" 
    
    setype = setype.lower();
           
    if len(sesize) != 2:
        raise StandardError('structureElement2D: sesize is not 2d');
        
    o = self.structureElementOffsets(sesize);
    omax = o.min(axis=1);
    sesize = numpy.array(sesize);
    
    if setype == 'sphere':
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1]];
        add = ((sesize + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];

        se = 1 - (x * x / (omax[0] * omax[0]) + y * y / (omax[1] * omax[1]));
        se[se < 0] = 0;
        return se / se.sum();

    elif setype == 'disk':
       
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1]];
        add = ((sesize + 1) % 2) / 2.;  
        
        x = g[0,:,:] + add[0];
        y = g[1,:,:] + add[1];

        se = 1 - (x * x / (omax[0] * omax[0]) + y * y / (omax[1] * omax[1]));
        se[se >= 0] = 1;
        se[se < 0] = 0;
        return se.astype('int');

    else:
        return numpy.ones(sesize);
    
    
    
def structureElement3D(setype = 'Disk', sesize = (3,3,3)):
    """Creates specific 3d structuring elements""" 
    
    setype = setype.lower();
           
    if len(sesize) != 3:
        raise StandardError('structureElement3D: sesize is not 3d');
            
    o = self.structureElementOffsets(sesize);
    omax = o.max(axis=1);
    sesize = numpy.array(sesize);
    
    if setype == 'sphere':
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1], -o[2,0]:o[2,1]];
        add = ((sesize + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];
        z = g[2,:,:,:] + add[2];

        se = 1 - (x * x / (omax[0] * omax[0]) + y * y / (omax[1] * omax[1]) + z * z / (omax[2] * omax[2]));
        se[se < 0] = 0;
        return se / se.sum();
        
    elif setype == 'disk':
        
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1], -o[2,0]:o[2,1]];
        add = ((sesize + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];
        z = g[2,:,:,:] + add[2];
        
        se = 1 - (x * x / (omax[0] * omax[0]) + y * y / (omax[1] * omax[1]) + z * z / (omax[2] * omax[2]));
        se[se < 0] = 0;
        se[se > 0] = 1;
        return se.astype('int');
    
    else:
        return numpy.ones(sesize);
        

def structureELement(setype = 'Disk', size = (3,3)):
    """Creates specific 2d and 3d structuring elements"""
    
    ndim = len(size); 
    if ndim == 2:
        return self.structureELement2D(setype, size);
    else:
        return self.structureElement3D(setype, size);
    
    
def filterKernel2D(ftype = 'Gaussian', size = (5,5), sigma = None, sigma2 = None, radius = None):
    """Creates a 2d filter kernel of a special type"""
    

    ftype = ftype.lower();
    o = self.structureElementOffsets(size);
    mo = o.min(axis=1);
    size = numpy.array(size);
    
    if ftype == 'mean':  # differnce of gaussians
        return numpy.ones(size)/ size.prod();
        
    elif ftype == 'gaussian':        
       
        if sigma == None:
           sigma = size / 2. / math.sqrt(2 * math.log(2));
           
        sigma = numpy.array(sigma);
        
        if len(sigma) < 3:
            sigma = numpy.array((sigma[0], sigma[0]));
        else:
            sigma = sigma[0:2];
                 
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1]];
        add = ((size + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];
        
        ker = numpy.exp(-(x * x / 2. / (sigma[0] * sigma[0]) + y * y / 2. / (sigma[1] * sigma[1])));
        return ker/ker.sum();
        
    elif ftype == 'sphere':
       
        if radius == None:
            radius = mo;
        radius = numpy.array(radius);
      
        if len(radius) < 3:
            radius = numpy.array((radius[0], radius[0]));
        else:
            radius = radius[0:2];
            
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1]];
        add = ((size + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];
        
        ker = 1 - (x * x / 2. / (radius[0] * radius[0]) + y * y / 2. / (radius[1] * radius[1]));
        ker[ker < 0] = 0.;
        return ker / ker.sum();
        
    elif ftype == 'disk':
       
        if radius == None:
            radius = mo;
        radius = numpy.array(radius);
        
        if len(radius) < 3:
            radius = numpy.array((radius[0], radius[0]));
        else:
            radius = radius[0:2];
            
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1]];
        add = ((size + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];
        
        ker = 1 - (x * x / 2. / (radius[0] * radius[0]) + y * y / 2. / (radius[1] * radius[1]));
        ker[ker < 0] = 0.;
        ker[ker > 0] = 1.0;
        return ker / ker.sum();
    
     
    elif ftype == 'log':  # laplacian of gaussians
       
        if sigma == None:
            sigma = size / 4. / math.sqrt(2 * math.log(2));
            
        sigma = numpy.array(sigma);

        if len(sigma) < 3:
            sigma = numpy.array((sigma[0], sigma[0]));
        else:
            sigma = sigma[0:2];

        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1]];
        add = ((size + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];

        ker = numpy.exp(-(x * x / 2. / (radius[0] * radius[0]) + y * y / 2. / (radius[1] * radius[1])));
        ker /= ker.sum();
        arg = x * x / math.pow(sigma[0], 4) + y * y/ math.pow(sigma[1],4) - (1/(sigma[0] * sigma[0]) + 1/(sigma[1] * sigma[1]));
        ker = ker * arg;
        return ker - ker.sum()/len(ker);

             
    elif ftype == 'dog':

        if sigma2 == None:
            sigma2 = size / 2. / math.sqrt(2 * math.log(2));
        sigma2 = numpy.array(sigma2);
        if len(sigma2) < 3:
            sigma2 = numpy.array((sigma2[0], sigma2[0]));
        else:
            sigma2 = sigma2[0:2];
        
        if sigma == None:
             sigma = sigma2 / 1.5;
        sigma = numpy.array(sigma);
        if len(sigma) < 3:
            sigma = numpy.array((sigma[0], sigma[0]));
        else:
            sigma = sigma[0:2];         
         
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1]];
        add = ((size + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];
        
        ker = numpy.exp(-(x * x / 2. / (sigma[0] * sigma[0]) + y * y / 2. / (sigma[1] * sigma[1])));
        ker /= ker.sum();
        sub = numpy.exp(-(x * x / 2. / (sigma2[0] * sigma2[0]) + y * y / 2. / (sigma2[1] * sigma2[1])));
        return ker - sub / sub.sum();
        
    else:
        raise StandardError('filter type ' + ftype + ' not implemented!');


def filterKernel3D(ftype = 'Gaussian', size = (5,5), sigma = None, sigma2 = None, radius = None):
    """Creates a 3d filter kernel of a special type"""
    
    ftype = ftype.lower();
    o = self.structureElementOffsets(size);
    mo = o.min(axis=1);
    size = numpy.array(size);
    
    if ftype == 'mean':  # differnce of gaussians
        return numpy.ones(size)/ size.prod();
        
    elif ftype == 'gaussian':        
       
        if sigma == None:
           sigma = size / 2. / math.sqrt(2 * math.log(2));
           
        sigma = numpy.array(sigma);
        
        if len(sigma) < 3:
            sigma = numpy.array((sigma[0], sigma[0], sigma[0]));
        else:
            sigma = sigma[0:3];
                 
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1], -o[2,0]:o[2,1]];
        add = ((size + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];
        z = g[2,:,:,:] + add[2];
        
        ker = numpy.exp(-(x * x / 2. / (sigma[0] * sigma[0]) + y * y / 2. / (sigma[1] * sigma[1]) + z * z / 2. / (sigma[2] * sigma[2])));
        return ker/ker.sum();
        
    elif ftype == 'sphere':
       
        if radius == None:
            radius = mo;
        radius = numpy.array(radius);
      
        if len(radius) < 3:
            radius = numpy.array((radius[0], radius[0], radius[0]));
        else:
            radius = radius[0:3];
            
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1], -o[2,0]:o[2,1]];
        add = ((size + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];
        z = g[2,:,:,:] + add[2];
        
        ker = 1 - (x * x / 2. / (radius[0] * radius[0]) + y * y / 2. / (radius[1] * radius[1]) + z * z / 2. / (radius[2] * radius[2]));
        ker[ker < 0] = 0.;
        return ker / ker.sum();
        
    elif ftype == 'disk':
       
        if radius == None:
            radius = mo;
        radius = numpy.array(radius);
      
        if len(radius) < 3:
            radius = numpy.array((radius[0], radius[0], radius[0]));
        else:
            radius = radius[0:3];
            
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1], -o[2,0]:o[2,1]];
        add = ((size + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];
        z = g[2,:,:,:] + add[2];
        
        ker = 1 - (x * x / 2. / (radius[0] * radius[0]) + y * y / 2. / (radius[1] * radius[1]) + z * z / 2. / (radius[2] * radius[2]));
        ker[ker < 0] = 0.;
        ker[ker > 0] = 1.0;
        return ker / ker.sum();
    
     
    elif ftype == 'log':  # laplacian of gaussians
       
        if sigma == None:
            sigma = size / 4. / math.sqrt(2 * math.log(2));
            
        sigma = numpy.array(sigma);

        if len(sigma) < 3:
            sigma = numpy.array((sigma[0], sigma[0], sigma[0]));
        else:
            sigma = sigma[0:3];

        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1], -o[2,0]:o[2,1]];
        add = ((size + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];
        z = g[2,:,:,:] + add[2];
        
        ker = numpy.exp(-(x * x / 2. / (radius[0] * radius[0]) + y * y / 2. / (radius[1] * radius[1]) + z * z / 2. / (radius[2] * radius[2])));
        ker /= ker.sum();
        arg = x * x / math.pow(sigma[0], 4) + y * y/ math.pow(sigma[1],4) + z * z / math.pow(sigma[2],4) - (1/(sigma[0] * sigma[0]) + 1/(sigma[1] * sigma[1]) + 1 / (sigma[2] * sigma[2]));
        ker = ker * arg;
        return ker - ker.sum()/len(ker);

             
    elif ftype == 'dog':

        if sigma2 == None:
            sigma2 = size / 2. / math.sqrt(2 * math.log(2));
        sigma2 = numpy.array(sigma2);
        if len(sigma2) < 3:
            sigma2 = numpy.array((sigma2[0], sigma2[0], sigma2[0]));
        else:
            sigma2 = sigma2[0:3];
        
        if sigma == None:
             sigma = sigma2 / 1.5;
        sigma = numpy.array(sigma);
        if len(sigma) < 3:
            sigma = numpy.array((sigma[0], sigma[0], sigma[0]));
        else:
            sigma = sigma[0:3];         
         
        g = numpy.mgrid[-o[0,0]:o[0,1], -o[1,0]:o[1,1], -o[2,0]:o[2,1]];
        add = ((size + 1) % 2) / 2.;
        x = g[0,:,:,:] + add[0];
        y = g[1,:,:,:] + add[1];
        z = g[2,:,:,:] + add[2];
        
        ker = numpy.exp(-(x * x / 2. / (sigma[0] * sigma[0]) + y * y / 2. / (sigma[1] * sigma[1]) + z * z / 2. / (sigma[2] * sigma[2])));
        ker /= ker.sum();
        sub = numpy.exp(-(x * x / 2. / (sigma2[0] * sigma2[0]) + y * y / 2. / (sigma2[1] * sigma2[1]) + z * z / 2. / (sigma2[2] * sigma2[2])));
        return ker - sub / sub.sum();
        
    else:
        raise StandardError('filter type ' + ftype + ' not implemented!');
        
        
def filterKernel(ftype = 'Gaussian', size = (5,5), sigma = None, radius = None, sigma2 = None):
    """Creates a filter kernel of a special type"""

    ndim = len(size);
    if ndim == 2:
        return self.filterKernel2D(ftype = ftype, size = size, sigma = sigma, sigma2 = sigma2, radius = radius);
    else:
        return self.filterKernel3D(ftype = ftype, size = size, sigma = sigma, sigma2 = sigma2, radius = radius);     
     
     
def hMaxTransform(img, h):
    #morphological reconstruction
    return reconstruction(img - h, img);
    
    
def regionalMax(img):
    return regmin(-img);
    
    
def extendedMax(img, h):
    #h max transform
    img = self.hMaxTransform(img, h);
        
    #regional max
    return self.regionalMax(img);
            
    
def detectCells(img, mask = None):
    """Detect Cells in 3d grayscale img using DoG filtering and maxima dtection"""
   
    # normalize data
    img = img.astype('float');
    dmax = 0.075 * 65535;
    ids = img > dmax;
    img[ids] = dmax;
    img /= dmax;    
  
    # background subtraction in each slice
    for z in range(img.shape[2]):
        img[:,:,z] = img[:,:,z] - morph.grey_opening(img[:,:,z], structure = self.structureELement('Disk', (30,30)));
        #img[:,:,z] = img[:,:,z] - morph.grey_opening(img[:,:,z], structure = self.structureELement('Disk', (150,150)));
        
    # mask
    if mask == None:
        mask = img > 0.01;
        mask = morph.binary_opening(mask, self.structureELement('Disk', (3,3,3)));
        

    #DoG filter
    fdog = self.filterKernel(ftype = 'DoG', size = [7,7,5]);
    imgdog = correlate(img, fdog);
    imgdog[imgdog < 0] = 0;
    imgdog /= imgdog.max();
    
    #iplt.plotTiling(imgdog) 

    imghmax = self.hMaxTransform(imgdog, 0.01);
    
    #iplt.plotTiling(imghmax) 
    
    imgmax = self.regionalMax(imghmax);
    
    #iplt.plotOverlayLabel(imghmax, imgmax.astype('int64'), alpha = False) 
    #iplt.plotOverlayLabel(img, imgmax.astype('int64'), alpha = True)     
    
    #extended maxima = h transform + regoinal maxima
    imgmax = self.extendedMax(imgdog, 0.01);
    imgmax = imgmax * mask;
 
    #iplt.plotOverlayLabel(img, imgmax.astype('int64'), alpha = False) 
    #iplt.plotOverlayLabel(imgdog, imgmax.astype('int64'), alpha = True) 

 
    #center of maxima
    imglab, nlab = label(imgmax);
    centers = numpy.array(center_of_mass(img, imglab, index = numpy.arange(1, nlab)));    

    
    #return centers, imglab, mask
    return centers;
    
    
from multiprocessing import Pool
#from functools import partial

def log_template(level, message):
    print("{}: {}".format(level, message))

 


#define the subroutine for the processing
def processSubStack(dsr):
    
    fn = dsr[0];
    segmentation = dsr[1];
    zrange  = dsr[2];
    iid = dsr[3];
    n = dsr[4];
    
    print "processing chunk " + str(iid) + "/" + str(n);
    
    f = io.openFile(fn, mode = 'r');
    dataset = io.readData(f, resolution=0);    
    
    #img = dataset[:,:, zrange[0]:zrange[1]];
    
    img = dataset[1200:1400,1200:1400, zrange[0]:zrange[1]];
    
    f.close();    
    
    return segmentation(img);
    
    
def parallelProcessStack(fn, chunksizemax = 100, chunksizemin = 30, chunkoverlap = 15, processes = 2, segmentation = self.detectCells):
    """Parallel segmetation on a stack"""
    
    f = io.openFile(fn);
    dataset = io.readData(f, resolution=0);
    
    #determine z ranges
    nz = dataset.shape[2];    
    nz = 200;    
    
    f.close();
    
    chunksize = chunksizemax;
    nchunks = int(math.ceil((nz - chunksize) / (1. * (chunksize - chunkoverlap)) + 1));   
    chunksize = (nz + (nchunks-1) * chunkoverlap) / nchunks;
    
    #increase overlap if chunks to small
    chunksizemin = min(chunksizemin, chunkoverlap);
    if chunksize < chunksizemin:
        print "parallelProcessStack: optimal chunk size " + str(chunksize) + " smaller than minimum chunk size " + str(chunksizemin) + "!"; 
        chunksize = chunksizemin;
        chunkoverlap = math.ceil(chunksize - (nz - chunksize) / (nchunks -1));
        print "parallelProcessStack: setting chunk overlap to " + str(chunkoverlap) + "!";
        
    #calucalte actual chunk sizes
    chunksizerest = chunksize;
    chunksize = math.floor(chunksize);
    chunksizerest = chunksizerest - chunksize;

    zranges = [(0, chunksize)];
    zcenters = [0];
    n = 1;
    csr = chunksizerest;
    zhi = chunksize;

    while (n < nchunks):
        n += 1;

        zhiold = zhi;
        zlo = zhi - chunkoverlap;
        zhi = zlo + chunksize;
        
        csr += chunksizerest;
        if csr > 1:
            csr = csr - 1;
            zhi += 1;
        
        zhi = min(zhi, nz);
        
        zranges.append((zlo, zhi));
        zcenters.append((zhiold - zlo) / 2. + zlo); 

    zcenters.append(nz);
    
    #print zcenters
    #print zranges    
    
    
    argdata = []
    for i in range(nchunks):
        argdata.append((f,segmentation, zranges[i], i, nchunks));
        
    
    # process in parallel
    pool = Pool(processes = processes);
    #print processes
    
    results = pool.map(processSubStack, argdata);
    
    #join the results
    for i in range(nchunks):
        cts = results[i];
        cts[:,2] += zranges[i][0];
        cts = cts[numpy.logical_and(cts[:,2] < zcenters[i+1], cts[:,2] >= zcenters[i]),:];
        results[i] = cts;
    
    return numpy.concatenate(results);

    
"""



# run segmentation
img = data[:,:,0:30];

import iDISCO.Segmentation.SpotDetection
import sys
self = sys.modules['iDISCO.Segmentation.SpotDetection'];



function hdImage = hdTransform2(I, h, n)
    disp('Enter into h-d transform');
    smallNumber = 1e-7;
    se = strel(ones(n));
    J = I-h;
    flag = 0;
%      [R C] = size(I);
%    count = 0;
    pI = J*0;
    while flag==0
        tempJ = min( imdilate(J, se), I);
        resud = sum( sum(abs(J-tempJ) ));
        J =  tempJ;
        pI = max(pI,J);
%        count = count +1
        if resud < smallNumber
            flag = 1;
        end
    end

    hdImage = I-pI;
%         figure
%     subplot(2,2,1)
%     imshow(I,[])
%     subplot(2,2,2)
%      imshow(pI,[])
%         subplot(2,2,3)
%      imshow( hdImage,[])
%      colormap('jet')
%      linkaxes;
    disp('Exit from h-d transform');
end 



    
    
import numpy as np
from skimage.morphology import reconstruction

x = np.linspace(0, 4 * np.pi)
y_mask = np.cos(x)

y_seed = y_mask.min() * np.ones_like(x)
y_seed[0] = 0.5
y_seed[-1] = 0
y_rec = reconstruction(y_seed, y_mask)
"""

    
    
    

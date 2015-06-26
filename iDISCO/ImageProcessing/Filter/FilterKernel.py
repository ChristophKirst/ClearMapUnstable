# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

"""
Filter Kernel

Routines to generate various filter kernels


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import numpy
import math

from iDISCO.ImageProcessing.Filter.StructureElement import structureElementOffsets


def filterKernel(ftype = 'Gaussian', size = (5,5), sigma = None, radius = None, sigma2 = None):
    """Creates a filter kernel of a special type"""
    
    ndim = len(size);
    if ndim == 2:
        return self.filterKernel2D(ftype = ftype, size = size, sigma = sigma, sigma2 = sigma2, radius = radius);
    else:
        return self.filterKernel3D(ftype = ftype, size = size, sigma = sigma, sigma2 = sigma2, radius = radius);


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
    o = structureElementOffsets(size);
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



def test():
    
        
    fk = self.filterKernel(ftype = 'Gaussian', size = (5,5), sigma = None, radius = None, sigma2 = None)
    
    
    
if __name__ == "__main__":
    
    self.test();
    
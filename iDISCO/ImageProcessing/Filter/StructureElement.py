# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

"""
Structure Elements

Routines to generate various structure elements


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import numpy


def structureElementOffsets(sesize):
    """Calculates offsets for a structural element given its size."""
    
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


def structureElement2D(setype = 'Disk', sesize = (3,3)):
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
        

def structureElement(setype = 'Disk', size = (3,3)):
    """Creates specific 2d and 3d structuring elements"""
    
    ndim = len(size); 
    if ndim == 2:
        return self.structureElement2D(setype, size);
    else:
        return self.structureElement3D(setype, size);
    
     
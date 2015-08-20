# -*- coding: utf-8 -*-
"""
iDISCO Settings 

Notes:
    - edit the indicated region to set up the ilastik and elastix paths

Created on Fri Aug 14 14:06:37 2015

@author: ckirst
"""

import sys
self = sys.modules[__name__];

import socket
import os

##############################################################################
# Path Settings 
##############################################################################

#path to ilastik 0.5 installation
IlastikPath = '/usr/local/ilastik-05-rc-final';

#path to eastix installation
ElastixPath = '/usr/local/elastix';


def setup():
    """Setup iDISCO for specific hosts"""
    global IlastikPath, ElastixPath
    
    hostname = socket.gethostname();
    
    if hostname == 'kagalaska.nld':  #Christophs Laptop 
        IlastikPath = '/home/ckirst/programs/ilastik-05/';
        ElastixPath = '/home/ckirst/programs/elastix/';
    
    elif hostname == 'mtllab-Ubuntu': #Nico workstation
        IlastikPath = '/usr/local/ilastik-05-rc-final';
        ElastixPath = '/usr/local/elastix';       
    
    ## insert your hostname specific settings here ##
    #elif hostname == 'your-host-name':
    #    IlastikPath = 'path-to-ilastik';
    #    ElastixPath = 'path-to-elastix';   
    ##

    # check existence:
    if not ElastixPath is None:
        if not os.path.exists(ElastixPath):
            raise RuntimeWarning('Settings: elastix path %s does not exists, cf. Settings.py');
            ElastixPath = None;
    
    if not IlastikPath is None:
        if not os.path.exists(IlastikPath):
            raise RuntimeWarning('Settings: ilastik path %s does not exists, cf. Settings.py');
            IlastikPath = None;

setup();


#path to iDISCO (used for testing and exmaples)
def iDISCOPath():
    """Returns path of IDISCO software"""
    fn = os.path.split(self.__file__)
    fn = os.path.abspath(fn[0]);
    return fn;

IDISCOPath = iDISCOPath();
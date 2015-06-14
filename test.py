# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 18:30:46 2015

@author: ckirst
"""

# -*- coding: utf-8 -*-
"""
Processing Brain Samples - Main Example Script


Created on Thu Jun  4 14:37:06 2015

@author: ckirst
"""

import numpy as np
import matplotlib as mpl

import iDISCO.IO.Imaris as io

from iDISCO.ImageProcessing.SpotDetection import detectCells
from iDISCO.ImageProcessing.ParallelProcessing import parallelProcessStack
from iDISCO.Visualization.Plot import plotTiling, plotOverlayLabel

from iDISCO.Utils.Timer import Timer


@profile
def run():
    timer = Timer();
    # open data
    fn = '/home/ckirst/Science/Projects/BrainActivityMap/Data/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims';
    #fn = '/run/media/ckirst/ChristophsBackuk4TB/iDISCO_2015_06/Adult cfos C row 20HF 150524.ims';
    #fn = '/home/nicolas/Windows/Nico/cfosRegistrations/Adult cfos C row 20HF 150524 - Copy.ims';
    #fn = '/home/ckirst/Science/Projects/BrainActivityMap/iDISCO_2015_04/test for spots added spot.ims'
    
    centers, centint = parallelProcessStack(fn, chunksizemax = 40, chunksizemin = 30, chunkoverlap = 15, 
                                            processes = 5, segmentation = detectCells, zrange = (800, 860));
                                            
    timer.printElapsedTime("Main");
    
    

if __name__ == '__main__':
    run()
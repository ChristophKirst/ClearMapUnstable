# -*- coding: utf-8 -*-
"""
*ClearMap* Registration, Image Analysis and Statistics Library.

*ClearMap* is a python toolbox for the analysis and registration of volumetric 
data from cleared tissues.

*ClearMap* is targeted towards large lightsheet volumetric imaging data
of iDISCO+ cleared mouse brains samples, their registration to the Alan brain atlas,
volumetric image processing and statistical analysis.

Author
""""""
     Christoph Kirst, 
     The Rockefeller University, New York City, 2015
   
License
""""""
     GNU GENERAL PUBLIC LICENSE Version 3

"""

__title__ = 'ClearMap'

__version__ = '0.9.2'

__author__ = 'Christoph Kirst'

__license__ = 'GNU GENERAL PUBLIC LICENSE Version 3'

__copyright__ = '2015 Christoph Kirst'

__all__ = ["Settings", "Parameter", "IO", "ImageProcessing", "Analysis"]

import ClearMap.Settings
import ClearMap.Parameter
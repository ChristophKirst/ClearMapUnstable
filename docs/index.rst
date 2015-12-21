ClearMap
========

*ClearMap* is a toolbox for the analysis and registration of volumetric data
from cleared tissues.

*ClearMap* has been designed to analyse large 3D image stack datasets obtained with Light Sheet Microscopy 
of iDISCO+ cleared mouse brains samples immunolabeled for nuclear proteins. ClearMap performs image registration to the Allen Institute Brain Atlas, volumetric image processing and object detection and statistical analysis. In particular, the tools in *ClearMap* have been written with the mapping of Immediate Early Genes in the brain as the primary application.

However, the tools included with *ClearMap* should also be more broadly useful for data obtained with other types of microscopes, other types of markers, and other clearing techniques. Moreover, the registration and region segmentation capabilities of ClearMap are not depending on the Atlases and annotations we used in our study. Users are free to import their own reference files and annotation files, so the use of *ClearMap* can be easily expanded to other species, other organs or samples.

*ClearMap* is written in Python 2.7, and is designed to take advantage of parallel processing capabilities of modern workstations. We hope the open structure of the code will enable in the future many new modules to be added to ClearMap to broaden the range of applications to different types of biological objects or structures.

Author and License
------------------

Authors:
^^^^^^^^

ClearMap lead programming and design:
"""""""""""""""""""""""""""""""""""""
Christoph Kirst, 
*The Rockefeller University*

Scripts and specific applications:
""""""""""""""""""""""""""""""""""
Nicolas Renier and Christoph Kirst
*The Rockefeller University*

Documentation:
""""""""""""""
Christoph Kirst and Nicolas Renier
*The Rockefeller University*

Additional programming and consulting:
""""""""""""""""""""""""""""""""""""""
Kannan Umadevi Venkataraju
*Cold Spring Harbor Laboratories*


License
^^^^^^^
GNU GENERAL PUBLIC LICENSE Version 3

See :download:`LICENSE <../LICENSE.txt>` or `gnu.org <http://www.gnu.org/licenses/gpl-3.0.en.html>`_ for details.




Using ClearMap
--------------

.. toctree::
   :maxdepth: 2

   introduction
   installation
   tutorial
   imageanalysis

ClearMap functions
------------------
.. toctree::
   :maxdepth: 2

   api/ClearMap 

Indices
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`






Introduction
============

*ClearMap* is a toolbox to analyze and register microscopy images of cleared 
tissue. It is targeted towards cleared brain tissue using the `iDISCO+ Clearing Method`_
but can be used with any volumetric imaging data. 

The package consists of four main modules:

   * `IO`_ for reading and writing images and data
   * `Alignment`_ for resampling, reorientation and registration of images onto references
   * `Image Processing`_ for correcting and quantifying the image data
   * `Analysis`_ for the statistical analysis of the data
   
IO
--

ClearMap supports a wide range of image formats with focus on volumetric data:

=========== ==========================================================
Format      Descrition
=========== ==========================================================
TIF         tif images and stacks
RAW / MHD   raw image files with optional mhd header file
NRRD        nearly raw raster data files
IMS         imaris image files
pattern     folder, file list or file pattern of a stack of 2d images
=========== ==========================================================

Images are represented internally as numpy arrays. CleaMap assumes images
in arrays are aranged as [x,y], [x,y,z] or [x,y,z,c] where x,y,z correspond to 
the x,y,z coordinates as when viewed in an image viewer such as [ImageJ]_ and
c to a possible color channel.

ClearMap also supports several data formats storing data points, such as
cell center coordinates or intensities

========= ==========================================================
Format    Descrition
========= ==========================================================
CSV       comma separated values in text file
NPY       numpy binary file
VTK       vtk point data file
========= ==========================================================

Alignment
---------

The Alignment module provides tools to resample, reorient and register
volumetric images in a fast parallel way.

Image registration is done by interfacing to the [Elastix]_ software package.

This package allows to align cleared mouse brains onto the Allan brain atlas [ABA]_.


Image Processing
----------------

ClearMap provides a number of image processing tools with focus on
processing large 3d volumetric images in parallel.


Main processing modules include:

======================================================= ===========================================================
Module                                                  Descrition
======================================================= ===========================================================
:mod:`~ClearMap.ImageProcessing.BackgroundRemoval`      Background estimation and removal via morphological opening
:mod:`~ClearMap.ImageProcessing.IlluminationCorrection` Correction of vignetting and other illumination errors
:mod:`~ClearMap.ImageProcessing.Filter`                 Filtering of the imagee via large set of filter kernels
:mod:`~ClearMap.ImageProcessing.GreyReconstruction`     Reconstruction of images
:mod:`~ClearMap.ImageProcessing.SpotDetection`          Detection of local peaks
:mod:`~ClearMap.ImageProcessing.CellDetection`          Detection of cell centers
:mod:`~ClearMap.ImageProcessing.CellSizeDetection`      Detection of cell shapes via watershed
:mod:`~ClearMap.ImageProcessing.IlastikClassification`  Classification of voxels via interface to [Ilastik]_ 
======================================================= ===========================================================

The modular structure of this sub-packages allows for fast and felxible integration of
additional modules.


Analysis
--------

This part of ClearMap provides a toolbox for the statistical analysis and 
visualization of detected cells or structures and region specific analysis
of annoated data.

For cleared mouse brains aligned to the [ABA]_ a wide range of statistical 
analysis tools with respect to the anotated brain regions in the atlas is
supported.

Key moduls are

====================================== =============================================================
Module                                 Descrition
====================================== =============================================================
:mod:`~ClearMap.Analysis.Voxelization` Voxelization of cells for visualization and analysis
:mod:`~ClearMap.Analysis.Statistics`   Statistical tools for the analysis of detected cells
:mod:`~ClearMap.Analysis.Label`        Tools to analysise data with espect to annotated refereneces
====================================== =============================================================


Examples
--------

Examples can be found :doc:`here </examples>` and in the api documentation.


iDISCO+ Clearing Method
-----------------------

iDISCO stands for immunolabeling-enabled three-dimensional imaging
of solvent-cleared organs. The iDISCO+ clearance method is 
decribed in [Renier2014]_ and [Renier2015]_.

iDISCO is modeled on classical histology techniques, facilitating translation 
of section staining assays to intact tissues, as evidenced by compatibility 
with 28 antibodies to both endogenous antigens and transgenic reporters 
like GFP. iDISCO enables facile volume imaging of immunolabeled structures 
in complex tissues

See these videos
   * `Dopaminergic system in the embryonic mouse <https://www.youtube.com/watch?v=-ctRUMQjizgvbLtLYkW6hI>`_
   * `Cortical and hippocampal neurons in the adult mouse brain <https://www.youtube.com/watch?v=vbLtLYkW6hI>`_

More info can be found on the [iDISCO]_ webpage.


References
----------

.. [Renier2014] `iDISCO: A Simple, Rapid Method to Immunolabel Large Tissue Samples
   for Volume Imaging, N. Renier, et al. 2014
   <http://dx.doi.org/10.1016/j.cell.2014.10.010>`_

.. [Renier2015] 'Mapping brain activity in the mouse at cellular resolution 
   with volume imaging using immediate early genes, N. Renier, et al. in prep.

.. [iDISCO] `iDISCO webpage, http://idisco.info/ < http://idisco.info/>`_

.. [ABA] `Allan Brain Atlas, http://www.brain-map.org/. <http://www.brain-map.org/>`_

.. [Elastix] `Elastix toolbox for rigid and nonrigid registration of 
   images, http://elastix.isi.uu.nl <http://elastix.isi.uu.nl>`_

.. [Ilastik] `Ilastik the interactive learning and segmentation toolkit, 
   http://ilastik.org/ <http://ilastik.org/>`_

.. [ImageJ] `ImageJ <http://imagej.net/Welcome>`_



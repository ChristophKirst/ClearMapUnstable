Installation
============

Requirements
^^^^^^^^^^^^

ClearMap is written in Python 2.7. It should run on any Python environment, but it also relies on external softwares such as `Elastix <http://elastix.isi.uu.nl>`_ which may not run optimally on Windows or Apple systems.
For typical use, we recommend a workstation running Ubuntu 14 or later with at least 4 CPU cores, 64Gb of RAM and SSD disks. 128Gb of RAM and 6 cores or above will have much increased performances. The processing time however will depend greatly on the parameters set, so your experience may be different. Also, large hard drives may be needed to host the raw data, although 1 to 4Tb of storage space should be enough for most users.

Installation
^^^^^^^^^^^^

To install ClearMap, first create a folder to contain all the files for the program. Then, open a terminal window, change to the directory hosting ClearMap and download the code from GIT by running this command:

    >>> git clone https://git.assembla.com/idisco.git

If you’re starting from scratch, you also need to download individually the following softwares:

  * To do the alignement, you should download `Elastix <http://elastix.isi.uu.nl>`_ (http://elastix.isi.uu.nl)

  * If you wish to use the machine learning filters, download `Ilastik <http://ilastik.org>`_ (http://ilastik.org), version 0.5 (this is the version implemented in ClearMap). This is an optional download, only if you wish to use this more complete object detection framework for complex objects.

And then the following libraries (most of them available from the software center of Ubuntu, or via the ``pip`` command or ``apt-get`` command in a terminal):

- Dev tools for Python 2.7 (from the software center)
- Matplotlib (from the software center)
- Numpy (from the software center)
- Scipy (from the software center)
- Skimage (from the software center)
- Mahotas (from the developer’s `website <http://luispedro.org/software/mahotas/>`_)
- h5py (from the software center) (for Imaris files input/output only)
- openCV (from the software center)
- MayaVI (from the software center)
- libboost-all-dev (from the software center)
- PyOpenGL (from the software center)
- qimage2ndarray (from the software center)
- PyQt4 (from the software center)
- tifffile (from the software center)
- EVTK (from the developer `website <https://bitbucket.org/pauloh/pyevtk/overview>`_) (only necessary, for output to vtk files)
- libhdf5-dev (from the software center)
- Cython (from the software center)

If you’re planning to use Ilastik, download `Vigra <http://ukoethe.github.io/vigra/>`_ (from the developer website).

We use `Spyder <https://pythonhosted.org/spyder/>`_ to run the code.

Configuration
^^^^^^^^^^^^^

Open the file ClearMap/Settings.py to set the paths of installations for Ilastik and Elastix:

    >>> IlastikPath = '/usr/local/ilastik-05-rc-final';
    >>> ElastixPath = '/usr/local/elastix';

Note that Ilastik is optional. If you haven’t installed it, you can set the path to ``None``. You can also set the installation to run on multiple machines by setting a host specific path:

    >>> if hostname == 'kagalaska.nld':  #Christoph’s Laptop 
    >>>         IlastikPath = '/home/ckirst/programs/ilastik-05/';
    >>>         ElastixPath = '/home/ckirst/programs/elastix/';
    >>>     elif hostname == 'mtllab-Ubuntu': #Nico’s Workstation
    >>>         IlastikPath = '/usr/local/ilastik-05-rc-final';
    >>>         ElastixPath = '/usr/local/elastix';       


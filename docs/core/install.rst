.. _Sec:install:

.. image:: images/SORA_logo.png
  :width: 500
  :align: center
  :alt: SORA: Stellar Occultation Reduction and Analysis

|
|


Installation
============


Python package requirements
---------------------------

Here we present the packages that will be used as base for our coding.
Those packages are also installed on the .

-  **astropy 4.0:** For astronomical related functions, mainly coordinates
   and time.

-  **astroquery 0.4.1:** To query astronomical database as JPL and Vizier.

-  **matplotlib 3.1.1:** For easy and beautiful plots.

-  **numpy 1.18.1:** Otimized mathematical functions.

-  **scipy 1.4.1:** Otimized functions for mathematics, science, and
   engineering.

-  **spiceypy 3.0.2:** SPICE/NAIF functions in python.

-  **pyerfa 2.0** Python wrapper for the ERFA library based on the SOFA library.   

-  **cartopy 0.17:** Geospatial data processing to produce maps.


Installing SORA
---------------

The user can install SORA and most of its requirements using **pip**, only
Cartopy should be installed by hand afterwards.

>>> user@path$> pip install sora-astro
>>> user@path$> conda install -c conda-forge cartopy

If you are a |GitHub| user, you can also use:

>>> user@path$> git clone https://github.com/riogroup/SORA/sora.git
>>> user@path$> cd sora
>>> user@path$> pip install .
>>> user@path$> conda install -c conda-forge cartopy

When new versions are available, the user can update it downloading the
last release from the SORA package in the riogroup organisation on
|GitHub|. If you want to be notified just follow the package.

.. |GitHub| raw:: html

   <a href="https://github.com/riogroup/SORA" target="_blank"> GitHub</a>


Functionalities
---------------

With SORA the user can:

#. Predict stellar occultations and obtain predictions maps;
#. Check when the stellar occultation will happen for a given observer;
#. Determine the immersion and emersion times from an lightcurve;
#. Check the chords in the skyplane;
#. Fit a circle for events with less than 3 chords and determine the 
   astrometrical position of the occulting object;
#. Fit an ellipse for events with more than 3 chords and determine the
   apparent size and shape of the occulting body and its position;

**All these steps can be found in our Jupyter-Notebooks Tutorials.**


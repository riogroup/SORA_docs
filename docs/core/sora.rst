.. _Sec:install:

.. image:: images/SORA_logo.png
  :width: 500
  :align: center
  :alt: SORA: Stellar Occultation Reduction and Analysis

|
|


Modules
^^^^^^^

SORA is a Python-based, object-oriented package for optimal analysis of
stellar occultation data. The user can use this package to build pipelines 
to analyse their stellar occultation’s data. Here follows the details for
each object class: **Body**, **Ephem**, **Star**, **PredictionTable**, 
**Observer**, **LightCurve**, **Occultation** and **ChiSquare**.


Body
----

The Body Class
==============
.. autoclass:: sora.Body
   :members:
   :private-members:

PhysicalData Class 
==================
.. autoclass:: sora.body.PhysicalData
   :members:

Module `utils`
==============
.. automodule:: sora.body.utils
   :members:


Ephem
-----

The EphemPlanete Class
======================
.. autoclass:: sora.EphemPlanete
   :members:

The EphemHorizons Class
=======================
.. autoclass:: sora.EphemHorizons
   :members:

The EphemKernel Class
=====================
.. autoclass:: sora.EphemKernel
   :members:

Module `utils`
==============
.. automodule:: sora.ephem.utils
   :members: 


Extra
-----

ChiSquare Class
===============
.. autoclass:: sora.extra.ChiSquare
   :members:

Plot ellipse
============
.. automodule:: sora.extra.plots
   :members:


Module `utils`
==============
.. automodule:: sora.extra.utils
   :members:


LightCurve
----------

The LightCurve Class
====================
.. autoclass:: sora.LightCurve
   :members:

Occutation Detection
====================
.. automodule:: sora.lightcurve.occdetect
   :members: occ_detect

Module's `utils`
================
.. automodule:: sora.lightcurve.utils
   :members: 


Observer
--------

The Observer Class
==================
.. autoclass:: sora.Observer
   :members:

Module `utils`
==============
.. automodule:: sora.observer.utils
   :members:


Occultation
-----------
The Occultation Class
=====================
.. autoclass:: sora.Occultation
   :members:

The Chord Class
===============
.. autoclass:: sora.occultation.Chord
   :members:

The Chordlist Class
===================
.. autoclass:: sora.occultation.chordlist.ChordList
   :members:

Fit Ellipse
===========
.. autofunction:: sora.occultation.fit_ellipse

Module `utils`
==============
.. automodule:: sora.occultation.utils
   :members:


Prediction
----------

Prediction Functions
====================
.. automodule:: sora.prediction.core
   :members:

Plot Occultation Map
====================
.. automodule:: sora.prediction.occmap
   :members:

PredictionTable Cbject
======================
.. autoclass:: sora.prediction.PredictionTable
   :members:


Star
----

The Star Class
==============
.. autoclass:: sora.Star
   :members:

Module `utils`
==============
.. automodule:: sora.star.utils
   :members:

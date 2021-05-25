.. _Sec:releases:

.. image:: images/SORA_logo.png
  :width: 500
  :align: center
  :alt: SORA: Stellar Occultation Reduction and Analysis

|
|


Releases
========

SORA v0.2 (Unreleased)
----------------------

.. warning::
    SORA v.02 is an unreleased version. The updates may contain errors that will be thoroughly tested 
    to assure the proper functionalities for each Class, Functions and Attributes.


New Features
^^^^^^^^^^^^

sora.body
~~~~~~~~~

- Created new Body Class which downloads the occulting body information from online source.
  At the moment, it downloads only from the Small-Body DataBase. The Body class will be the manager
  for all the Body informations, such as Ephem, Shape, Ring, etc. [:issue:`51`]

- New Class PhysicalData, which inherits from astropy.units.quantity.Quantity, is created to handle
  physical data with uncertainty, reference and notes. [:issue:`51`]

- "pole_position_angle" and "apparent_magnitude" functions are now present in Body instead 
  of Ephem. [:issue:`51`]

- Created a hardcoded satellite database to complement missing data of SBDB. It must be replaced 
  in the future. [:issue:`61`]

sora.config
~~~~~~~~~~~

sora.ephem
~~~~~~~~~~

- A new EphemHorizons was created which is strictly equal to EphemJPL (EphemJPL may be removed 
  in v1.0). [:issue:`51`]

sora.extra
~~~~~~~~~~

- Allow two ChiSquare objects to be combined into one: `chi3 = chi1 + chi2`. [:issue:`61`]

- New function get_ellipse_points() that calculates the positions on the perimeter of an 
  ellipse. [:issue:`60`]

sora.lightcurve
~~~~~~~~~~~~~~~

sora.observer
~~~~~~~~~~~~~

sora.occultation
~~~~~~~~~~~~~~~~

- A shortcut was created in Occultation where the user can pass the coordinate of the star 
  directly to Occultation, the Star object will be created automaticaly. [:issue:`46`]

- New Chord Class introduced to handle a chord with an Observer and a LightCurve. [:issue:`53`]

- New ChordList Class introduced to handle the list of Chords in an Occultation. [:issue:`53`]

- New function .get_impact_param() that calculatesthe impact parameter, minimal distance
  between the chord and the centre position, in Chord and ChordList. [:issue:`60`]

- New function .get_theoretical_times(),that calcultates the theoretical times and chord size
  for a given ellipse in Chord and ChordList. [:issue:`60`]

- New fucntion .check_time_shift() that calculates the offset in time to align the center of 
  the chords in Occultation. [:issue:`60`]

- New parameters sigma_result, that saves the result with an extended error bar, and 
  ellipse_error, that adds a further systematic error to be considered, in 
  Occultation.fit_ellipse(). [:issue:`60`]

- New function fiter_negative_chord() that compares the ChiSquare from an Ellipse fitting 
  with the chords and remove the solutions that would cross a negative chord [:issue:`60`]

sora.prediction
~~~~~~~~~~~~~~~

- prediction() now makes use of the user input of the star to calculate faster the 
  occultation parameters. [:issue:`48`]

- prediction() now can make predictions using Gaia-EDR3. A new parameter "catalogue" was created
  for choosing between Gaia-DR2 and Gaia-EDR3. [:issue:`61`]

- Fixed bug when plotting the heights in the map in a rotated projection. [:issue:`54`]

sora.star
~~~~~~~~~


API Changes
^^^^^^^^^^^

- Update the argument "log" to "verbose" on all modules. [:issue:`61`]

sora.config
~~~~~~~~~~~

sora.ephem
~~~~~~~~~~

- "pole_position_angle" and "apparent_magnitude" is passed to Body Class. In Ephem, it will raise
  a FutureWarning. [:issue:`51`]

- The Ephem classes are now passed through the Body Class which will have priority over Ephem
  attributes. Parameters such as "spkid", "radius", "H" and "G". [:issue:`51`]

- All Ephem Classes now inherits from BaseEphem, which holds core functionality for all of 
  them. [:issue:`51`]

sora.extra
~~~~~~~~~~

sora.lightcurve
~~~~~~~~~~~~~~~

- Removed the necessity for LightCurve to have a unique name associated. [:issue:`53`]

- Cycle time is now determined via mode instead of median. [:issue:`56`]

sora.observer
~~~~~~~~~~~~~

- Removed the necessity for Observer to have a unique name associated. [:issue:`53`]

sora.occultation
~~~~~~~~~~~~~~~~

- The new Body Class was implemented in Occultation. For backward compatibility, the previous
  usage is still possible if the Ephem object have a name. The Body Class is only required
  if the object is a planet or a planetary satellite. [:issue:`51`]

- Deprecated some functions that were passed to ChordList. [:issue:`53`]

sora.prediction
~~~~~~~~~~~~~~~

- prediction() now creates the time array inside each division to avoid memory overflow. 
  [:issue:`48`]

- prediction() now propagates the positions of the stars using only the proper motions
  before comparing the stars with the ephemeris. [:issue:`48`]

- The new Body Class was implemented in prediction. For backward compatibility, the previous
  usage is still possible. [:issue:`51`]

sora.star
~~~~~~~~~


Bug Fixes
^^^^^^^^^

sora.config
~~~~~~~~~~~

sora.ephem
~~~~~~~~~~

sora.extra
~~~~~~~~~~

sora.lightcurve
~~~~~~~~~~~~~~~

- Corrected bug in LightCurve model where the size of the star was being interpreted
  as radius instead of diameter. [:issue:`60`]

sora.observer
~~~~~~~~~~~~~

sora.occultation
~~~~~~~~~~~~~~~~

sora.prediction
~~~~~~~~~~~~~~~

- Fixes issue that happenned in occ_params() when the instant of the occultation was outside the 
  given range. The function now gives appropriate error messages. The automatic range search was 
  increased to 50 min from central instant in a recursive search. [:issue:`45, 48`]

sora.star
~~~~~~~~~


SORA v0.1.2 (2020/Dec/14)
-------------------------

New Features
^^^^^^^^^^^^

sora.config
~~~~~~~~~~~

sora.ephem
~~~~~~~~~~

sora.extra
~~~~~~~~~~

sora.lightcurve
~~~~~~~~~~~~~~~

sora.observer
~~~~~~~~~~~~~

sora.occultation
~~~~~~~~~~~~~~~~

sora.prediction
~~~~~~~~~~~~~~~

sora.star
~~~~~~~~~

- Star() is now able to fully receive astrometric parameters from the user. [:issue:`48`]

- Star() is able to download and use the distance from Bailer-Jones et al (2018). [:issue:`27`]

- Gaia-EDR3 was implemented in Star() and is now a default feature. [:issue:`52`]


API Changes
^^^^^^^^^^^

sora.config
~~~~~~~~~~~

sora.ephem
~~~~~~~~~~

sora.extra
~~~~~~~~~~

sora.lightcurve
~~~~~~~~~~~~~~~

sora.observer
~~~~~~~~~~~~~

sora.occultation
~~~~~~~~~~~~~~~~

sora.prediction
~~~~~~~~~~~~~~~

sora.star
~~~~~~~~~

- The star module was moved to its own directory. [:issue:`52`]


Bug Fixes
^^^^^^^^^

sora.config
~~~~~~~~~~~

sora.ephem
~~~~~~~~~~

sora.extra
~~~~~~~~~~

sora.lightcurve
~~~~~~~~~~~~~~~

sora.observer
~~~~~~~~~~~~~

sora.occultation
~~~~~~~~~~~~~~~~

sora.prediction
~~~~~~~~~~~~~~~

sora.star
~~~~~~~~~

- Star now calculates the robust propagation of the position of the star and 
  correspondent uncertainties. [:issue:`18`]

- Fixed bug in Star().__str__() where pmDEC was printed wrong. [:issue:`43`]

- A small bug fix was made in Star with the units of the star position error 
  when coordinates are local. [:issue:`51`]


SORA v0.1.1 (2020/Jul/30)
-------------------------

New Features
^^^^^^^^^^^^

sora.config
~~~~~~~~~~~

- Module to verify if kwargs are allowed was created. This was included throughout the code. 
  [:issue:`8`]

sora.ephem
~~~~~~~~~~

sora.extra
~~~~~~~~~~

- Added a parameter that allows the used to plot a dot corresponding
  the center of the ellipse [:issue:`35`]

sora.lightcurve
~~~~~~~~~~~~~~~

- Property LightCurve.time_mean that returns the mean time of the chord (positive) or
  the mean time of the observation (negative). [:issue:`34`]

sora.observer
~~~~~~~~~~~~~

- Function Observer.altaz() that calculates the altitude and azimuth for a given target 
  and instant. [:issue:`34`]

sora.occultation
~~~~~~~~~~~~~~~~

sora.prediction
~~~~~~~~~~~~~~~

- Four new parameters were added to `plot_occ_map()`: `path`: for the user to select
  a directory where to save the plots; `site_name`: If True, the name of the sites
  will be plotted; `chord_delta` and `chord_geo`: for the user to plot the path of
  a chord from distance of the center or passing by some coordinate, respectively. [:issue:`35`]

- Two methods were added to `PredictionTable()` to help the user to remove bad events
  from table: `keep_from_selected_images()` and `remove_occ()`. [:issue:`35`]

sora.star
~~~~~~~~~


API Changes
^^^^^^^^^^^

sora.config
~~~~~~~~~~~

- config module is now a directory. It now includes a module with decorators
  and another for variables. [:issue:`31, 35`]

sora.ephem
~~~~~~~~~~

- In EphemKernel, `code` argument was replaced by `spkid`. When using 'code',
  a FutureWarning is raised stating `code` as deprecated and will be removed from v1.0. [:issue:`26`]

sora.extra
~~~~~~~~~~

sora.lightcurve
~~~~~~~~~~~~~~~

- In LightCurve.immersion and LightCurve.emersion, an error will rise when these values were not 
  instanciated or fitted. [:issue:`34`]

- Now the user has the possibility to redefine `tref`, `immersion`, `emersion`,
  `initial_time` and `end_time` after instantiated. [:issue:`35`]

- `lambda_0` argument was replaced by `central_bandpass` and `delta_lambda` by `delta_bandpass`. 
  When using 'lambda_0' or `delta_lambda`, a FutureWarning is raised stating `lambda_0` or 
  `delta_lambda` as deprecated and will be removed from v1.0. [:issue:`36`]

sora.observer
~~~~~~~~~~~~~

sora.occultation
~~~~~~~~~~~~~~~~

- Occultation.new_astrometric_positions() now shows a warning when time is far
  by more than 1 day from the occultation closest approach. [:issue:`21`]

- Occultation.to_log() and print(Occultation) added the polar radius, equivalent radius, 
  the Sun-Geocenter-Target angle and the Moon-Geocenter-Target angle, geocentric albedo,
  the altitude and azimuth of the target for each Observer. [:issue:`17`]

- In `fit_ellipse()`, `pos_angle` and `dpos_angle` were deprecated in favor of
  `position_angle` and `dposition_angle`. [:issue:`35`]

- Changed "GCRS" to "Geocentric" in the string representation to avoid confusion
  about the reference frame. [:issue:`35`]
  
sora.prediction
~~~~~~~~~~~~~~~

- prediction() now calculates the ephemeris inside each division to avoid memory overflow. 
  [:issue:`31`]

- PredictionTable.to_ow() will now raise a warning if the radius or the error of
  the ephemeris is not present. [:issue:`35`]

sora.star
~~~~~~~~~

- Now Star downloads all parameters from Gaia and saves them in the `meta_gaia` attribute 
  [:issue:`35`]


Bug Fixes
^^^^^^^^^

sora.config
~~~~~~~~~~~

sora.ephem
~~~~~~~~~~

- Added function get_position() to EphemPlanete. This corrects a bug that prevented
  Occultation to run with EphemPlanete. [:issue:`41`]

- Fixed bug in EphemJPL where `id_type` was redefined inside __init__(). [:issue:`41`]

sora.extra
~~~~~~~~~~

sora.lightcurve
~~~~~~~~~~~~~~~

- Fixed error that appears when the fit was done separately (immersion and emersion times). 
  Now the final model agrees with the fitted values. [:issue:`9`]

- Fixed error when the file with the light curve has three columns. [:issue:`19`]

- Fixed error when the exptime within the LightCurve was set as zero or negative. [:issue:`23`]

- Fixed error in the automatic mode of LightCurve.normalize(). [:issue:`34`]

- Fixed bug that was raised in LightCurve.log() when there were no initial or end times
  for lightcurves instantiated with immersion and emersion. [:issue:`35`]

sora.observer
~~~~~~~~~~~~~

sora.occultation
~~~~~~~~~~~~~~~~

- Corrected error calculation using err = sqrt(star_err^2 + fit_err^2) [:issue:`18`]

- Occultation.plot_occ_map() now uses the fitted ellipse to calculate the projected shadow 
  radius [:issue:`22`]

- Corrected bug that raised an error when calling Occultation.get_map_sites()
  and there were no observation added to Occultation. [:issue:`31`]

- Corrected bug that did not save the fitted params in all occultations when
  more than one occultation was used in fit_ellipse(). [:issue:`35`]

- Added `axis_labels` and `lw` (linewidth) to Occultation.plot_chords(). [:issue:`35`]

sora.prediction
~~~~~~~~~~~~~~~

- Fixed error that was generated when only one prediction was found. [:issue:`16`]

- Fixed error in the output format of PredictionTable.to_ow() when coordinate was positive 
  [:issue:`35`]

sora.star
~~~~~~~~~


SORA v0.1 (2020/May/20)
-----------------------

Object Classes
^^^^^^^^^^^^^^

The documentation of all classes and functions are on their docstrings, while the 
scientific part is presented in the full documentation. Here follows a list with the 
main Object Classes:

Ephem 
~~~~~

Three Object Classes created to generate geocentric ephemeris for a given solar 
system object. **EphemJPL** queries the |JPL Horizons| service and download ephemeris 
information. **EphemKernel** reads the BSP files to calculate the ephemeris using the Spiceypy 
package. **EphemPlanet** reads an ASCII file with previously determined positions and interpolate 
them for a given instant.

.. |JPL Horizons| raw:: html

   <a href="https://ssd.jpl.nasa.gov/horizons.cgi" target="_blank">JPL Horizons</a>


Star 
~~~~

Object Class created to deal with the star parameters. 
From the Gaia-DR3 or DR2 Source ID or a sky region, it queries the VizieR 
service and downloads the star’s information. From Gaia Catalog 
(Gaia Collaboration et al., |Gaia2016a|, |Gaia2016b|, |Gaia2018|, |Gaia2020|) it 
gets the RA, DEC, parallax, proper motions, G magnitude and star radius; 
from the NOMAD Catalog (Zacharias et al. |NZ04|) it gets the B, V, R, J, H and 
K magnitudes. The user can calculate the ICRS coordinate of the star at any epoch. 
It can be barycentric (corrected from proper motion) or geocentric (corrected from 
proper motion and parallax). Also, the apparent diameter of the star is calculated 
using Gaia information, or some models such as Van Belle (|VB99|) and  
Kervella et al. (|PK04|).

.. |Gaia2016a| raw:: html

   <a href="https://ui.adsabs.harvard.edu/abs/2016A\%26A...595A...1G/abstract" target="_blank">2016a</a>

.. |Gaia2016b| raw:: html

   <a href="https://ui.adsabs.harvard.edu/abs/2016A\%26A...595A...2G/abstract" target="_blank">2016b</a>

.. |Gaia2018| raw:: html

   <a href="https://ui.adsabs.harvard.edu/abs/2018A\%26A...616A...1G/abstract" target="_blank">2018</a>

.. |Gaia2020| raw:: html

   <a href="https://ui.adsabs.harvard.edu/abs/2020arXiv201201533G" target="_blank">2020</a>

.. |VizieR| raw:: html

   <a href="https://vizier.u-strasbg.fr/viz-bin/VizieR" target="_blank">VizieR</a>

.. |NZ04| raw:: html

   <a href="https://ui.adsabs.harvard.edu/abs/2004AAS...205.4815Z/abstract" target="_blank">2004</a>

.. |VB99| raw:: html

   <a href="https://ui.adsabs.harvard.edu/abs/1999PASP..111.1515V/abstract" target="_blank">1999</a>

.. |PK04| raw:: html

   <a href="https://ui.adsabs.harvard.edu/abs/2004A%26A...426..297K/abstract" target="_blank">2004</a>


Observer
~~~~~~~~
Object Class created to deal with the observer location. The user can also download 
the ground-based observatories from the Minor Planet Center (|MPC|) database.

.. |MPC| raw:: html

   <a href="https://minorplanetcenter.net/iau/lists/ObsCodesF.html" target="_blank">MPC</a>

Light Curve
~~~~~~~~~~~

Object Class that receives the observational light curve (with time and the occulted 
star normalized photometry relative to reference stars) and some observational 
parameters (filter and exposure time). It has functions to determine the instants 
that the solar system object enters in front of the star and leaves, (immersion and 
emersion times, respectively). The model considers a sharp-edge occultation model 
(geometric) convolved with Fresnel diffraction, stellar diameter (projected at the body 
distance) and finite integration time (Widemann et al., |Wi09|; Sicardy et al., |BS11|).

.. |Wi09| raw:: html

   <a href="https://ui.adsabs.harvard.edu/abs/2009Icar..199..458W/abstract" target="_blank">2009</a>

.. |BS11| raw:: html

   <a href="https://ui.adsabs.harvard.edu/abs/2011Natur.478..493S/abstract" target="_blank">2011</a>

Occultation
~~~~~~~~~~~

Main Object Class within SORA, created to analyze stellar occultations, and control all 
the other Object Classes within this package. Its functions allow converting the times 
for each observatory in the occulted body positions in the sky plane relative to the 
occulted star (f, g) (|IERS|). Also, to obtain the best ellipse parameters (centre 
position, apparent equatorial radius, oblateness and the position angle of the apparent 
polar radius) that fit the points. The results are the apparent size, shape and astrometrical 
position of the occulting body.

.. |IERS| raw:: html

   <a href="https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html" target="_blank">IERS Conventions</a>


Some extra Objects Classes:

PredictionTable
~~~~~~~~~~~~~~~

Using the **prediction** function within SORA results in an Object Class that is a 
slight modification of an AstropyTable. The added changes allow to create the occultation 
map for each prediction, convert into specific formats, such as |OccultWatcher| and 
PRAIA (Assafin et al. (|MA11|)).

.. |OccultWatcher| raw:: html

   <a href="https://www.occultwatcher.net/" target="_blank">OccultWatcher</a>

.. |MA11| raw:: html

   <a href="https://ui.adsabs.harvard.edu/abs/2011gfun.conf...85A/abstract" target="_blank">2011</a>


ChiSquare
~~~~~~~~~

This Object Class is the result of the fitting functions within SORA, such as 
*LightCurve.occ_lcfit()* and *Occultation.fit_ellipse()*. This Class has functions that 
allow viewing the values that minimize the :math:`\chi^2` tests, the uncertainties within 
:math:`n\sigma`, plotting the tests, and saving the values.   


Input and outputs
^^^^^^^^^^^^^^^^^

Inputs
~~~~~~

  - **Event Related (Star and Ephem)**
 
    - Object Name or provisory designation
    - Object Code (only for EphemKernel)
    - BSP file and name (only for EphemKernel)
    - DE file and name (only for EphemKernel)
    - Ephemeris offset for RA and DEC - :math:`\Delta \alpha \cdot \cos \delta`,
      :math:`\Delta \delta` (set as 0,0)
    - Occultation date and time
    - Occulted star coordinates RA and DEC; or Gaia code
    - Star offset for RA and DEC - :math:`\Delta \alpha \cdot \cos \delta`, 
      :math:`\Delta \delta` (set as 0,0)

  - **Observer Related**
 
    - Site name and location (latitude, longitude, and height; or IAU/MPC code)
    - Light curve file and name; or array with fluxes and times; or immersion 
      and emersion times
    - Exposure time in seconds
    - Observational bandwidth in microns (set as 0.7 +/- 0.3 microns, Clear)

  - **Fitting Related**
 
    - Initial guess for light curve fitting: immersion, emersion and opacity.
    - Range to explore all three parameters
    - Initial guess for ellipse parameters: center (f,g), equatorial radius, 
      oblateness, and position angle
    - Range to explore all five parameters


Outputs
~~~~~~~

  - **Star**
 
    - Star Gaia-DR2 ID
    - Star coordinates at 2015.5 and uncertainty - RA and DEC (hh mm ss.sss , 
      +dd mm ss.sss, mas, mas)
    - Star proper motion - in RA, DEC - and uncertainties (mas/yr)
    - Star parallax and uncertainty (mas)
    - Star coordinates propagated to event epoch and uncertainty - RA and DEC 
      (hh mm ss.sss , +dd mm ss.sss, mas, mas)
    - Star magnitudes G, B, V, R, J, H, K (mag)
    - Star projected diameter and model (km and mas, model: GDR2, Van Belle, Kervella)
    - Star offset applied in RA and DEC (mas, mas)


  - **Object and Ephemeris**

    - Object Name
    - Object radius (km)
    - Object mass (kg)
    - Ephemeris kernel (version and DE)
    - Offset applied in RA/DEC (mas, mas)
    - Object’s distance (AU)
    - Object apparent magnitude for the date (mag)

  - **Occultation**

    - Event date and time (yyyy-mm-dd hh:mm:ss.sss)
    - Closest approach Angle - CA (arcsec)
    - Reference time (yyyy-mm-dd hh:mm:ss.sss)
    - Position Angle - PA (degree)
    - Shadow’s velocity relative to the geocenter (km/s)
    - Number of positive observations
    - Number of negative observations


  - **Observer Information**
 
    - Detection status (positive, negative, overcast, tech. problem, other)
    - Site Name
    - Site MPC/IAU code (if any)
    - Site coordinates - Latitude, Longitude and height  (dd mm ss.s ; dd mm ss.s ; m)
    - Light curve file name
    - Number of images (lines in LC)

  - **Light curve fitting information (for each positive detection)**

    - Acquisition start time (yyyy-mm-dd hh:mm:ss.sss)
    - Acquisition end time (yyyy-mm-dd hh:mm:ss.sss)
    - Exposure time (s)
    - Cycle time (s)
    - Time offset applied in LC (s)
    - Light curve calculated RMS
    - Calculated normalised flux and bottom flux (standard = 1, 0)
    - Band width and uncertainty (microns)
    - Shadow's velocity relative to the station (km/s)
    - Fresnel scale (s and km)
    - Projected stellar size scale (s and km)
    - Integration time scale (s and km)
    - Dead time scale (s and km)
    - Model resolution - size of synthetic LC point (s and km)
    - Immersion Time and uncertainty (yyyy-mm-dd hh:mm:ss.sss +/- s.sss)
    - Immersion Time and uncertainty - :math:`1\sigma` and :math:`3\sigma` (s)
    - Emersion Time and uncertainty (yyyy-mm-dd hh:mm:ss.sss +/- s.sss)
    - :math:`\chi^2` fit model
    - Emersion Time and uncertainty - :math:`1\sigma` and :math:`3\sigma` (s)
    - Minimum Chi-square - :math:`\chi^2_{min}`
    - Number of fitted points for im- and emersion
    - Number of fitted parameters
    - Minimum Chi-square per degree of freedom - :math:`\chi^2_{min-pdf}`

  - **Elipse fit procedure**
 
    - Fitted parameters: Equatorial radius and uncertainty (km); Center position 
      :math:`(f_0, g_0)` and :math:`1\sigma` uncertainties (km, km); Oblateness and 
      uncertainty; Position angle and uncertainty (degree)
    - Minimum Chi-square -  :math:`\chi^2_{min}`
    - Minimum Chi-square per degree of freedom - :math:`\chi^2_{min-pdf}`
    - Number points used to fit ( X points from Y chords )
    - Astrometric object center position at occ. time and uncertainty (hh mm ss.sss 
      +dd mm ss.sss +/- mas)

  - **Plots and files (some are optional)**

    - Prediction map (Lucky Star model)
    - Normalised light curve - for each site :math:`(x = time; y = flux)`
    - Chi-square map for immersion and emersion times :math:`(x = time; y = \chi^2)`
    - Light curve and synthetic LC- for each site :math:`(x = time; y = flux)`
    - Chords projected in sky plane :math:`(x = \xi (km); y = \eta (km) )`
    - Chi-square map for each ellipse parameter :math:`(x = time; y = \chi^2_{param})`
    - Chords projected in sky plane and the best ellipse fitted with :math:`1\sigma` 
      uncertainties :math:`(x = \xi (km); y = \eta (km) )`
    - Log file with all information


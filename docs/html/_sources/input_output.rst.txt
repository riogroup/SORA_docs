.. _Sec:input_output:

Input and Outputs
=================

.. _SubSec_input_list:

INPUT list
----------

Here it is listed the inputs that the user needs to provide for SORA
v0.1 to work.

#. **Event Related (Star and Ephem)**

   -  Object Name or provisory designation

   -  Object Code / SPKID (only for EphemKernel)

   -  BSP file and name (only for EphemKernel)

   -  DE file and name (only for EphemKernel)

   -  Ephemeris offset for RA and DEC -
      :math:`\Delta \alpha \cdot \cos \delta`, :math:`\Delta \delta`
      (set as 0,0)

   -  Occultation date and time

   -  Occulted star coordinates RA and DEC; or Gaia code

   -  Star offset for RA and DEC -
      :math:`\Delta \alpha \cdot \cos \delta`, :math:`\Delta \delta`
      (set as 0,0)

#. **Observer Related**

   -  Site name and location (latitude, longitude, and height; or
      IAU/MPC code)

   -  Light curve file and name; or array with fluxes and times; or
      immersion and emersion times

   -  Exposure time in seconds

   -  Observational bandwidth in microns (set as 0.7 :math:`\pm` 0.3
      microns, Clear)

#. **Fitting Related**

   -  Initial guess for light curve fitting: immersion, emersion and
      opacity.

   -  Range to explore all three parameters

   -  Initial guess for ellipse parameters: center (f,g), equatorial
      radius, oblateness, and position angle

   -  Range to explore all five parameters

OUTPUT list
-----------

Here we present all the information that is obtained through
calculations, from queries, or from plots that are performed in SORA
code.

#. **Generic Information** (for a specific occultation):

   -  **Star**

      -  Star Gaia-DR2 ID

      -  Star coordinates at 2015.5 and uncertainty - RA and DEC (hh mm
         ss.sss , +dd mm ss.sss, mas, mas)

      -  Star proper motion - in RA, DEC - and uncertainties (mas/yr)

      -  Star parallax and uncertainty (mas)

      -  Star coordinates propagated to event epoch and uncertainty - RA
         and DEC (hh mm ss.sss , +dd mm ss.sss, mas, mas)

      -  Star magnitudes G, B, V, R, J, H, K (mag)

      -  Star projected diameter and model (km and mas, model: GDR2, Van
         Belle, Kervella)

      -  Star offset applied in RA and DEC (mas, mas)

   -  **Object and Ephmeris**

      -  Object Name

      -  Object radius (km)

      -  Object mass (kg)

      -  Ephemeris kernel (version and DE)

      -  Offset applied in RA/DEC (mas, mas)

      -  Object distance (AU)

      -  Object apparent magnitude for the date (mag)

      -  Object absolute magnitude

      -  Object geometric albedo (V-band)

   -  **Occultation**

      -  Event date and time (yyyy-mm-dd hh:mm:ss.sss)

      -  Closest approach Angle - CA (arcsec)

      -  Reference time (yyyy-mm-dd hh:mm:ss.sss)

      -  Position Angle - PA (degree)

      -  Shadow’s velocity relative to the geocenter (km/s)

      -  Number of positive observations

      -  Number of negative observations

#. **Observer Information** (for each observer):

   -  Detection status (positive, negative, overcast, tech. problem,
      other)

   -  Site Name

   -  Site MPC/IAU code (if any)

   -  Site coordinates - Latitude, Longitude and height (dd mm ss.s ; dd
      mm ss.s ; m)

   -  Light curve file name

   -  Number of images (lines in LC)

#. **Light curve fitting information** (for each positive detection):

   -  Acquisition start time (yyyy-mm-dd hh:mm:ss.sss)

   -  Acquisition end time (yyyy-mm-dd hh:mm:ss.sss)

   -  Exposure time (s)

   -  Cycle time (s)

   -  Time offset applied in LC (s)

   -  Light curve calculated RMS

   -  Calculated normalised flux and bottom flux (standard = 1, 0)

   -  Band width and uncertainty (microns)

   -  Shadow’s velocity relative to the station (km/s)

   -  Fresnel scale (s and km)

   -  Projected stellar size scale (s and km)

   -  Integration time scale (s and km)

   -  Dead time scale (s and km)

   -  Model resolution - size of synthetic LC point (s and km)

   -  Immersion Time (yyyy-mm-dd hh:mm:ss.sss)

   -  Immersion Time uncertainty - 1\ :math:`\sigma` and
      3\ :math:`\sigma` (s)

   -  Emersion Time (yyyy-mm-dd hh:mm:ss.sss)

   -  :math:`\chi^2` fit model

   -  Emersion Time uncertainty - 1\ :math:`\sigma` and
      3\ :math:`\sigma` (s)

   -  Minimum Chi-square - :math:`\chi^2_{min}`

   -  Number of fitted points for im- and emersion

   -  Number of fitted parameters

   -  Minimum Chi-square per degree of freedom -
      :math:`\chi^2_{min-pdf}`

#. **Elipse fit procedure:**

   -  Fitted parameters: Equatorial radius and uncertainty (km); Center
      position (:math:`f_0`, :math:`g_0`) and 1\ :math:`\sigma`
      uncertainties (km, km); Oblateness and uncertainty; Position angle
      and uncertainty (degree)

   -  Minimum Chi-square - :math:`\chi^2_{min}`

   -  Minimum Chi-square per degree of freedom -
      :math:`\chi^2_{min-pdf}`

   -  Number points used to fit ( X points from Y chords )

   -  Astrometric object center position at occ. time and uncertainty
      (hh mm ss.sss +dd mm ss.sss :math:`\pm` mas)

#. **Plots and files:** (some are optional)

   -  Prediction map (Lucky Star model)

   -  Normalised light curve - for each site (x = time; y = flux)

   -  Chi-square map for immersion and emersion times (x = time; y =
      :math:`\chi^2`)

   -  Light curve and synthetic LC- for each site (x = time; y = flux)

   -  Chords projected in sky plane (x = :math:`\xi` (km); y =
      :math:`\eta` (km) )

   -  Chi-square map for each ellipse parameter (x = time; y =
      :math:`\chi^2_{param}`)

   -  Chords projected in sky plane and the best ellipse fitted with
      1\ :math:`\sigma` uncertainties (x = :math:`\xi` (km); y =
      :math:`\eta` (km) )

   -  Log file with all information
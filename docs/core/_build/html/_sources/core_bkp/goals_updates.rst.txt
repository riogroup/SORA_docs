.. _Sec:goals_updates:

Goals and Version updates
=========================

-  **Documentation** [SORAversion0.0]

   -  Program design

   -  Pipeline workflow

   -  Input and Output logs and Figures

   -  SORA description

   -  Jupyter notebook with a complete reduction.

   -  Jupyter notebook tutorial for Ephem().

   -  Jupyter notebook tutorial for Star().

   -  Jupyter notebook tutorial for Observer().

   -  Jupyter notebook tutorial for LightCurve().

   -  Jupyter notebook tutorial for Occultation().

   -  Jupyter notebook tutorial for Prediction().

   -  Jupyter notebook tutorial for ChiSquare().


-  **Version 0.1**

   -  Calculate ephemeris from BSP kernels with SPICE.

   -  Calculate ephemeris querying the JPL Horizons website.

   -  Calculate ephemeris (using :math:`\xi`, :math:`\eta`) in the same
      sense as ephem_planete,

   -  Calculate position angle of the pole.

   -  Download star coordinates, magnitudes, proper motions, parallax
      and uncertainties from Gaia-DR2.

   -  Download star magnitudes B, V, R, J, H, K and uncertainties from
      NOMAD.

   -  Propagate star coordinates to event epoch (proper motion and
      parallax).

   -  Propagate star coordinates uncertainties (first order
      approximation).

   -  Predict stellar occultations and calculate their parameters.

   -  Include ephemeris uncertainties for occultation prediction
      (instead of a fixed value of object radius + Earth Radius).

   -  Download object apparent magnitude from JPL.

   -  Calculate object apparent magnitude from absolute magnitude.

   -  Calculate Magnitude drop and expected bottom flux.

   -  Calculate star apparent size using equations from van Belle
      (1999).

   -  Calculate star apparent size using equations from Kervella (2004).

   -  Calculate star apparent size using Gaia diameter and object
      distance.

   -  Search for site(s) coordinates in the MPC database (geocentric).

   -  Read light curve file with 2 columns (time and flux).

   -  Light curve normalisation using a polynomial fit.

   -  Automatically detect stellar occultation(s) in the light curve
      (determination of initial guess for immersion and emersion times).

   -  Create a synthetic stellar occultation light curve.

   -  :math:`\chi^2` minimisation and fit of immersion and emersion
      times from light curve

   -  :math:`\chi^2` minimisation and fit of opacity from light curve

   -  Analyse chi-square and determine the values that minimise it.

   -  Calculate the projected points :math:`f_i` and :math:`g_i` in the
      plane of sky of each occultation instant :math:`i` (immersions and
      emersions) with uncertainties.

   -  Fit ellipses using the projected points (f, g).

   -  Calculate event velocity for the site with respect to the normal
      of the ellipse surface (current method assumes a circle).

   -  Calculate object RA and DEC for the ellipse fitted center.

   -  Complete output log file

   -  Output figures and plots

   -  Fully integrated SORA modules.

   -  Routines adapted for user input (offline version): star, ephem,
      site(s) coordinates (topocentric).


-  **Version 0.1.1**

   -  Bugs and Error fixes.

   -  Inputs and outputs list update.

   -  Documentation and “README” file update.

   -  PEP8 verification

   -  Warnings messages and a few Error bypass fixes.

   -  Automatic calculation of object geometric albedo from object
      absolute magnitude.

   -  Read light curve files with User format

   -  Include Solar and Lunar separation angle in log file.

   -  Include observation elevation and azimuth (in degrees) in log
      file.


-  **Version 0.X - Expected short and medium term implementations**
   [SORAversion0.2]

   -  Bugs and Error fixes. [Under analysis]

   -  Include gravitational light shift for shadow path. [Under
      analysis]

   -  Include uncertainty in star diameter from GAIA. [Under analysis]

   -  Automatic initial value for ellipse fit. [Under analysis]

   -  Download object rotation period. [Under analysis]

   -  Use rotation period to determine density for Maclaurin. [Under
      analysis]

   -  Use rotation period to determine density for Jacobi. [Under
      analysis]

   -  Develop a F-test to verify any indication of atmosphere. [Under
      analysis]

   -  Improvement for secondary event determination. [Under analysis]

   -  Routine to estimate lower limits for presence of atmosphere.
      [Future work]

   -  Query to database with sites of recurrent amateurs observers.
      [Future work]

   -  From pole position angle determine expected pole coordinates
      (ICRS). [Under development]

   -  Determine true longitude in the equatorial plane (ring plane) of a
      point in the sky plane (:math:`f_i,g_i`) [Under development]

   -  Determine radial velocity for a point in the sky plane
      (:math:`f_i,g_i`) in the equatorial plane (ring plane). [Under
      development]

   -  Routine to read SORA output logfile as an input for SORA (avoid
      recalculations). [Future work]

   -  GUI front end design. [Future work]

   -  Model for integration with LIneA TNO portal (user access control).
      [Future work]


-  **Version 1.X and future updates - Expected long term implementations**

   -  Bugs and Error fixes

   -  Fit a 3D model to the occultation data.

   -  Density estimation for differentiated objects (see Hely’s
      dissertation)

   -  Implementation of different :math:`\chi^2` fitting procedures
      (e.g. fit 2 or 4 degree polynom in :math:`\chi^2` values for
      immersion and emersion instants).

   -  Fresnell diffraction functions depending on pre-set values for
      camera sensitivity and filter bandpass.

   -  Create a “reverse engeneering” process to determine a possible
      detection in a observing site.

   -  Star diameter bin calculation (current: star size divided by 24,  
      see Sec. `Light Curve fit <#SubSec:LC_fit>`__).

   -  New reduction routines for objects with atmosphere

   -  Include module to calculate object rotation period,

   -  Integration with google maps for prediction maps (creating kmz
      files),

   -  Integration with the stellar occultation database,
      (http://occultations.ct.utfpr.edu.br/)

   -  SORA available for private use (core team)

   -  SORA available for public use (occultation collaborators)

   -  SORA available for public use (general public under conditions)

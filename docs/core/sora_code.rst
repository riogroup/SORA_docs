
.. _SubSec:code:

SORA Code
=========

.. _SubSec:code_packages:

Python package requirements
---------------------------

Here we present the packages that will be used as base for our coding.
Those packages are also installed on the .

-  astropy 4.0: For astronomical related functions, mainly coordinates
   and time.

-  astroquery 0.4.1: To query astronomical database as JPL and Vizier.

-  cartopy 0.17: Geospatial data processing to produce maps.

-  matplotlib 3.1.1: For easy and beautiful plots.

-  numpy 1.18.1: Otimized mathematical functions.

-  scipy 1.4.1: Otimized functions for mathematics, science, and
   engineering.

-  spiceypy 3.0.2: SPICE/NAIF functions in python.

.. _SubSec:code_package_install:

Python package installation
---------------------------

Note that on LIneA docker all the packages will be installed and the
user do not need to worry about packages installation. You only need to
install if you want to use SORA in your personal computer.

The best way to install Python 3.7 is using anaconda:
https://www.anaconda.com/products/individual. They have versions for
different Operational System (Windows, MacOS and Linux). After
installing Anaconda, to verify which packages and version you have
installed, open a terminal window (in Unix systems [7]_) and type:

::

       user@path$> conda list 

or to search for a specific package:

::

       user@path$> conda list | grep <package>

As an example, to search for the scipy package:

::

       user@path$> conda list | grep scipy

In order to install a specific Python package, in a terminal window
type:

::

       user@path$> conda install <package>

or to install a specific version of the package:

::

       user@path$> conda install <package>=<version>

As an example, to install the scipy package version 1.4.1:

::

       user@path$> conda install scipy=1.4.1

If conda can not install a package, the user can try using *pip* instead
of *conda*, although it is recommended always using the *conda*.

For the specific case of SORA, the following commands should be executed
for installing all requirements [8]_.

::

       user@path$>conda install astropy=4.0 cartopy=0.17 spiceypy=3.0 matplotlib scipy ipython jupyter
       user@path$> conda install -c conda-forge cartopy
       user@path$> pip install --pre astroquery

.. _SubSec:sora_install:

SORA installation
-----------------

After installing all the python packages (from the previous session), to
install SORA locally, you simply need to clone the git repository, and
then install the package:

::

       user@path$> git clone https://github.com/riogroup/SORA/sora.git
       user@path$> cd sora
       user@path$> pip install .

.. _SubSec:sora_update:

SORA Update
-----------

When new versions are available, the user can update it downloading the
last release from the SORA package in the riogroup organisation on
GitHub. If you want to be notified just follow the package.

.. _SubSec:code_classes:

Object Classes
--------------

Here we present the structured object classes that SORA contains. In the
subsections we provide Jupyter-Notebooks containing the main functions
within each Object Class.

-  **Ephem:** Three Object Classes created to generate geocentric
   ephemeris for a given solar system object. **EphemJPL** queries the
   JPL Horizons service and download ephemeris information.
   **EphemKernel** reads the BSP files to calculate the ephemeris using
   the Spiceypy package. **EphemPlanet** reads an ASCII file with
   previously determined positions and interpolate them for a given
   instant.

-  **Star:** Object Class created to deal with the star parameters. From
   the Gaia-DR2 Source ID or a sky region, it queries the VizieR service
   and downloads the star’s information. From Gaia DR2 Catalog (Gaia
   Collaboration 2016a, 2016b and 2018) it gets the RA, DEC, parallax,
   proper motions, G magnitude and star radius; from the NOMAD Catalog
   (Zacharias et al. 2004) it gets the B, V, R, J, H and K magnitudes.
   The user can calculate the ICRS coordinate of the star at any epoch.
   It can be barycentric (corrected from proper motion) or geocentric
   (corrected from proper motion and parallax). Also, the apparent
   diameter of the star is calculated using Gaia DR2 information, or
   some models such as Van Belle (1999) and Kervella et al. (2004).

-  **Observer:** Object Class created to deal with the observer
   location. The user can also download the ground-based observatories
   from the Minor Planet Center (MPC) database.

-  **Lightcurve:** Object Class that receives the observational light
   curve (with time and the occulted star normalized photometry relative
   to reference stars) and some observational parameters (filter and
   exposure time). It has functions to determine the instants that the
   solar system object enters in front of the star and leaves,
   (immersion and emersion times, respectively). The model considers a
   sharp-edge occultation model (geometric) convolved with Fresnel
   diffraction, stellar diameter (projected at the body distance) and
   finite integration time (Widemann et al., 2009; Sicardy et al.,
   2011).

-  **Occultation:** Main Object Class within SORA, created to analyze
   stellar occultations, and control all the other Object Classes within
   this package. Its functions allow converting the times for each
   observatory in the occulted body positions in the sky plane relative
   to the occulted star (:math:`f`, :math:`g`) (IERS Conventions). Also,
   to obtain the best ellipse parameters (centre position, apparent
   equatorial radius, oblateness and the position angle of the apparent
   polar radius) that fit the points. The results are the apparent size,
   shape and astrometrical position of the occulting body.

-  **PredictionTable:** Using the *prediction()* function within SORA
   results in an Object Class that is a slight modification of an
   AstropyTable. The added changes allow to create the occultation map
   for each prediction, convert into specific formats, such as
   OccultWatcher and PRAIA (Assafin et al., 2011).

-  **ChiSquare:** This Object Class is the result of the fitting
   functions within SORA, for instance the *LightCurve.occ_lcfit()* and
   *Occultation.fit_ellipse()*. This Class has functions that allow
   viewing the values that minimize the :math:`\chi^2` tests, the
   uncertainties within :math:`n-\sigma`, plotting the tests, and saving
   the values.

.. _SubSubSec:code_classes_extra:

Module: extra
~~~~~~~~~~~~~

-  Set of independent functions that are used by *Observer* but do not
   depend on it.

   -  **draw_ellipse()** – Plots an ellipse given input parameters.

      **INPUT**: *equatorial_radius*: Semi-major axis of the ellipse.
      *oblateness*: Oblateness of the ellipse. Default=0.0. *center_x*
      and *center_y*: Coordinate of the ellipse (x,y respectively).
      Default=0.0. *pos_angle*: Pole position angle. Default=0.0.
      *\**kwargs*: all other parameters will be parsed directly by
      matplotlib.

-  **ChiSquare** – Object to handle :math:`\chi^2`.

   -  **\__init__()** – Instantiate a ChiSquare object.

      **INPUT**: *chi2*: An array with with the values of chisquare.
      *npts*: Number of points used in the fit. This is not the number
      of attempts. *\**kwargs*: Every new argument will be interpreted
      as data. The name of the argument is the key to access values.
      They must have the same size as chi2. At least one kwargs must be
      given.

      .. code:: python

         from sora.extra import ChiSquare

         chi1 = ChiSquare(chi2=chi_arr, npts=10, immersion=immersion_arr, emersion=emersion_arr)
         chi2 = ChiSquare(chi2=chi_arr, npts=10, radius=many_values, oblateness=other_many_values)

   -  **get_nsigma()** – get mean values and range within given sigma.

      **INPUT**: *sigma*: sigma to calculate range of values; *key*: if
      no key is given, it calculates for every param. if given, it
      calculates only for given param.

      **OUTPUT**: if no key is given, it returns a dictionary with the
      valeus for all params. if key is given, it return an array with
      the mean value and error bar.

   -  **get_values()** – Similar to get_nsigma, but it returns arrays
      with all the values for each key

      **INPUT**: *sigma*: sigma to get range of values.

   -  **plot_chi2()** – plots param x chi2.

      **INPUT**: *key*: if key is given, it will only plot for this
      parameter, otherwise, it will plot for every param.

   -  **to_file()** – Saves data in a file

      **INPUT**: *namefile*: Name of the file to save data. It will also
      save another file named namefile.label stating what each column
      is.

.. _SubSubSec:code_classes_occobs:

Module: observer
~~~~~~~~~~~~~~~~

-  Set of independent functions that are used by *Observer* but do not
   depend on it.

   -  **search_code_mpc()** – Read page
      https://www.minorplanetcenter.net/iau/lists/ObsCodes.html which
      contains all observers in the MPC, which is frequently updated. It
      contains the IAU code, longitude, latitude, altitude, and name for
      the sites.

      **INPUT**: None

      **OUTPUT**: python dictionary with information of all available
      sites in the MPC database.

   .. code:: python

      from sora.observer import search_code_mpc

      observatories = search_code_mpc()
      #(out) Looking the MPC database ...

      print(observatories)
      """
      {'000': ('Greenwich', <EarthLocation (0.62411, 0., 0.77873) earthRad>),
       '001': ('Crowborough', <EarthLocation (0.62991772, 0.0016953, 0.77411) earthRad>),
       ...
       'Z99': ('Clixby Observatory, Cleethorpes', <EarthLocation (0.59546796, -0.00022095, 0.800687) earthRad>)}
       """

-  **Observer**

   -  **\__init__()** - instantiate the observer object.

      **INPUT**: IAU code to search for the sites in the MPC database or
      name, latitude, longitude and altitude;

      **OUTPUT** Site name, longitude, latitude and altitude.

      .. code:: python

         from sora import Observer

         obs = Observer(code='874')
         print(obs.name)
         # 'Observatorio do Pico dos Dias, Itajuba'

         from astropy.coordinates import EarthLocation
         site = EarthLocation('-45 34 57', ' -22 32 04', 1864)
         lna = Observer(name='OPD', site=site)

         lna = Observer(name='OPD', lon='-45 34 57', lat='-22 32 04', height=1864)

         print(lna)
         #Site: OPD
         #Geodetic coordinates: Lon: -45d34m57s, Lat: -22d32m04s, height: 1.864 km

   -  **sidereal_time()** - Calculates the sidereal time for the time
      and local given.

      **INPUT**: *time*: an Astropy.time object; *mode*: location to
      calculate the sidereal time that can be ’greenwich’ or ’local’
      (apparent for the site location – default).

      **OUTPUT**: sidereal time for the local chosen in *mode*.

      .. code:: python

         from sora import Observer
         from astropy.time import Time

         obs = Observer(code='874')
         obs.sidereal_time('2019-06-07 03:54:22.60', 'local')
         # 17h53m05.8251s
         obs.sidereal_time('2019-06-07 03:54:22.60', 'greenwich')
         # 20h55m25.6251s

   -  **get_ksi_eta()** - Calculates the :math:`\xi` and :math:`\eta`
      from input time and projected at the direction of star.

      **INPUT**: *time*: an Astropy.time object or string in the ISO
      format; *star*: The coordinate of the star in the same frame as
      the ephemeris. It can be a string in the format “hh mm ss.s +dd mm
      ss.ss” or an astropy SkyCoord object.

      **OUTPUT**: :math:`\xi`, :math:`\eta`.

      .. code:: python

         from sora import Observer

         casleo = Observer('Casleo', '-69 17 44.9', '-31 47 55.6', 2492)
         casleo.get_ksi_eta(time='2019-06-07 03:54:22.60', star='19 21 18.63201 -21 44 25.3924')
         # -3911.0928016639878 -1705.5362129689515

Disclosure: Once an Observer is set, no other Observer object with the
same name can be set again until the first is deleted; After an Observer
is set, the user can correct any value directly to the Observer
attribute, for example:

.. code:: python

   from sora import Observer

   lna = Observer(name='OPD', lon='45 34 57', lat='-22 32 04', height=1864)
   print(lna)
   #Site: OPD
   #Geodetic coordinates: Lon: 45d34m57s, Lat: -22d32m04s, height: 1.864 km

   lna.lon = '-45 34 57'
   print(lna)
   #Site: OPD
   #Geodetic coordinates: Lon: -45d34m57s, Lat: -22d32m04s, height: 1.864 km

.. _SubSubSec:code_classes_star:

Module: star
~~~~~~~~~~~~

-  Set of independent functions that are used by *Star* but do not
   depend on it.

   -  **search_star()** – Search for stars on the Vizier catalogues:

      **INPUT**: code of the star or coord and radius for search, the
      name of the columns to download, the catalogue from where to
      download.

      **OUTPUT**: Astropy table with all the stars.

      .. code:: python

         from sora.star import search_star
         import astropy.units as u
         columns = ['Source', 'RA_ICRS', 'e_RA_ICRS', 'DE_ICRS', 'e_DE_ICRS', 'Plx', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Gmag', 'e_Gmag', 'Dup', 'Epoch', 'Rad']
         cat = search_star(code="4117746607441803776", columns=columns, catalog='I/345/gaia2')
         cat2 = search_star(coord='19 21 18.63201 -21 44 25.3924', radius=2*u.arcsec, columns=columns, catalog='I/345/gaia2')
         print(cat2)
         #TableList with 1 tables:
         #   '0:I/345/gaia2' with 15 column(s) and 1 row(s) 
         #       Source           RA_ICRS     e_RA_ICRS ... Dup Epoch   Rad  
         #                          deg          mas    ...       yr    Rsun 
         #------------------- --------------- --------- ... --- ------ ------
         #6772694935064300416 290.32764199494    0.0385 ...   0 2015.5   1.10

   -  **van_belle()** – Calculates the star diameter in mas using the
      equations from van Belle (1999).

      **INPUT**: magB, magV, magK.

      **OUTPUT**: python dictionary with diameter for all magnitudes and
      star types.

   -  **kervella()** – Calculates the star diameter in mas using the
      equations from Kervella (2004).

      **INPUT**: magB, magV, magK.

      **OUTPUT**: python dictionary with diameter for all magnitudes.

      .. code:: python

         from sora.star import van_belle, kervella
         van_belle(magB=10, magV=11, magK=12)
         #{'sg': {'B': <Quantity 0.01614359 mas>, 'V': <Quantity 0.01761976 mas>},
         # 'ms': {'B': <Quantity 0.00831764 mas>, 'V': <Quantity 0.01086426 mas>},
         # 'vs': {'B': <Quantity 0.02254239 mas>, 'V': <Quantity 0.02685344 mas>}}
         kervella(magB=10, magV=11, magK=12)
         #{'V': <Quantity 0.01100272 mas>, 'B': <Quantity 0.01020704 mas>}

-  **Star**

   -  **\__init__()** – instantiate the star object.

      **INPUT**: gaia code (Source) or coordinate. It will update the
      coordinate from Gaia and read the magnitudes from NOMAD.

      .. code:: python

         from sora import Star

         star = Star(code='4117746607441803776') # Gaia-DR2 Source code
         star2 = Star(coord='19 21 18.63201 -21 44 25.3924')
         print(star2)
         #ICRS star coordinate at J2015.5: RA=19h21m18.63408s +/- 0.0385 mas, DEC=-21d44m25.3767s +/- #0.0364 mas
         #Gaia-DR2 star Source ID: 6772694935064300416
         #Magnitudes: G: 15.292, B: 14.380, V: 14.870, R: 15.150, J: 14.119, H: 13.783, K: 13.685
         #Diameter: 0.0081 mas, Source: Gaia-DR2

   -  **\__searchgaia()** – hidden function that searchs for the star in
      the Gaia catalogue and saves the information. It is called once in
      the \__init__(). If the star radius is not found it gives a
      warning.

   -  **\__getcolors()** – hidden function that searchs for the star in
      the NOMAD catalogue and saves the magnitudes. If some band is not
      found, it gives a warning.

   -  **set_magnitude()** – Set the magnitudes of a star in any band as
      given by the user. It saves the magnitude in a python dictionary.

      **INPUT**: any ‘band=value‘, see example below.

      .. code:: python

         star.set_magnitude(G=10)
         # UserWarning: G mag already defined. G=15.292 will be replaced by G=10
         star.set_magnitude(X=25.0)
         print(star.mag)
         #{'G': 10, 'X': 25.0}

   -  **set_diameter()** – Set the user diameter of the star in mas
      which has higher priority than other diameters.

   -  **van_belle()** – call van_belle function passing the star
      magnitudes.

   -  **kervella()** – call kervella function passing the star
      magnitudes.

   -  **apparent_diameter()** – calculate the apparent diameter of the
      star given a distance furnish by the user.

      **INPUT**: distance to calculate apparent radius (required); mode,
      the function to calculate, where it can be ‘*user*’, ‘*gaia*’,
      ‘*kervella*’, ‘*van_belle*’ or ‘*auto*’, where ‘*auto*’ will run
      all the functions in this order until it is able to obtain a
      value. Default mode=‘*auto*’. In case mode=‘*kervella*’, the user
      must give parameter obs_filter that can be ‘B’ or ‘V’. In case
      mode=‘*van_belle*’, the user must give parameters obs_type and
      star_type, that can be ‘*sg*’ (super giant), ‘*ms*’ (main
      sequence) or ‘*vs*’ (variable star).

      .. code:: python

         star.apparent_diameter(10)
         #Apparent diameter using Gaia
         #1.287568km

         star.apparent_diameter(5, mode='kervella', obs_filter='V')
         #Apparent diameter using Kervella et al. (2004)
         #0.86978164km

         star.apparent_diameter(5, mode='van_belle', obs_filter='V', star_type='ms')
         #Apparent diameter using van Belle (1999)
         #0.75472762km

   -  **barycentric()** – returns the position of the star corrected
      from proper motion in the ICRS.

      **INPUT**: *time*: Instant of observation.

      **OUTPUT**: Astropy SkyCoord Object in the barycentric reference
      corrected from proper motion.

   -  **geocentric()** – returns the position of the star corrected from
      proper motion and parallax in the GCRS.

      **INPUT**: *time*: Instant of observation.

      **OUTPUT**: Astropy SkyCoord Object in the geocentric reference
      corrected from proper motion and parallax.

      .. code:: python

         star.geocentric('2019-06-07 03:54:22.60').to_string('hmsdms', precision=5)
         # 19h21m18.63201s -21d44m25.39241s
         star.barycentric('2019-06-07 03:54:22.60').to_string('hmsdms', precision=5)
         # 19h21m18.63198s -21d44m25.39247s

.. _SubSubSec:code_classes_ephem:

Module: ephem
~~~~~~~~~~~~~

In the file ‘**ephem.py**’, there are four classes:
‘\ *EphemPlanete*\ ‘, ‘\ *EphemKernel*\ ‘, ‘\ *EphemJPL*\ ‘ and
‘\ *Ephemeris*\ ‘. The first, ‘EphemPlanete‘ emulates what’s done by the
ephem_planete and fit_d2_ksi_eta fortran codes. The second ‘EphemKernel‘
makes use of spiceypy package (python wrapper of SPICE C functions) to
calculate the ephemeris. ‘EphemJPL‘ queries the JPL Horizons website for
the ephemeris desired. The last one, ‘Ephemeris’ class would be used as
controller. If a text file is given, it points to EphemPlanete,
otherwise if it is given *bsp* files, it points to other classes, and so
on. For now, it only points to ‘EphemPlanete‘.

-  Set of independent functions that are used by *Ephemeris* objects but
   do not depend on it.

   -  **read_obj_data()** – reads data table from
      http://devel2.linea.gov.br/~altair.gomes/radius.txt with the
      radius and ephem error for selected objects.

-  **EphemPlanete** – Class that controls the ephemeris simulating the
   fortran program ephem_planete

   -  **\__init__()** – reads the file with the ephemeris, saves the
      coordinates and time.

      **INPUT**: the name of the object and the name of the file with
      the ephemeris. The file must have JD, RA and DEC in degrees, and
      object geocentric distance in AU.

   -  **fit_d2_ksi_eta()** – calculates the on-sky difference between
      star and ephemeris (:math:`\xi,\eta`) and fits a two degree
      function.

      **INPUT**: an Astropy SkyCoord object.

   -  **get_ksi_eta()** – returns the calculated :math:`\xi` and
      :math:`\eta` from input time, the star coordinate can be passed as
      well to skip the call of fit_d2_ksi_eta().

      **INPUT**: time (necessary), star coordinate (not necessary).

      **OUTPUT**: pair of ksi and eta. If the input time is an array,
      two arrays for ksi and eta will be returned.

Example of how to use EphemPlanete

.. code:: python

   from sora import EphemPlanete

   ephem = EphemPlanete('Phoebe', 'ephem_phoebe_ph15.txt')

   ephem.fit_d2_ksi_eta(star='19 21 18.63201 -21 44 25.3924')
   #Fitting ephemeris position relative to star coordinate 19h21m18.632s -21d44m25.3924s
   #ksi = aksi*t2 + bksi*t + cksi
   #eta = aeta*t2 + beta*t + ceta
   #t=(jd-2458641.62083333)/(2458641.70416667-2458641.62083333)
   #        aksi=74.33578049015848
   #        bksi=120036.66960513973
   #        cksi=-56728.3072442839
   #        aeta=11.57934951952548
   #        beta=17144.759007406785
   #        ceta=-6972.117976312244
   #Residual RMS: ksi=0.004 km, eta=0.001 km

   ephem.get_ksi_eta(time=Time('2019-06-07 03:54:22.60'))
   # 3685.627507205092 1657.0083431480398

   print(ephem)
   #Ephemeris of Phoebe.
   #Valid from 2019-06-07 02:54:00.000 until 2019-06-07 04:54:00.000

-  **EphemKernel** – Class that controls the ephemeris given BSP
   kernels.

   -  **\__init__()** – reads the BSP files given

      **INPUT**: *name*: the name of the object; *code*: the code in the
      kernel of the target; *kernels*: the list of kernels. The list of
      kernels must be in the order of SPICE reading, with the last one
      having highest priority. For instance, if the same path is present
      in two given kernels, the last one will be used.

   -  **get_position()** – Calculates and return the geocentric position
      of the target for the given time.

      **INPUT**: time (necessary).

      **OUTPUT**: Astropy SkyCoord Object with the geocentric coordinate
      of the object for the given times.

   -  **get_ksi_eta()** – returns the calculated :math:`\xi` and
      :math:`\eta` from input time, the star coordinate can be passed as
      well to skip the call of fit_d2_ksi_eta().

      **INPUT**: time (necessary), star coordinate (not necessary).

      **OUTPUT**: pair of :math:`\xi` and :math:`\eta`. If the input
      time is an array, two arrays for :math:`\xi` and :math:`\eta` will
      be returned.

   -  **get_pole_position_angle()** – Calculates and return the position
      and the aperture angle of the target’s pole for the given time.

      **INPUT**: pole: Pole coordinate (necessary); time: (necessary).

      **OUTPUT**: position angle of the pole in degrees; aperture angle
      in degrees.

.. code:: python

   from sora import EphemKernel
   from astropy.time import Time

   ephe2 = EphemKernel('Phoebe', 609, kernels=['path/ph15.bsp', 'path/de438.bsp'])

   ephe2.get_position('2019-06-07 03:54:22.60')

   ksi, eta = ephe2.get_ksi_eta(time='2019-06-07 03:54:22.60', star='19 21 18.63201 -21 44 25.3924')

   ephe2.get_pole_position_angle(pole='10 05 12.000 +41 28 48.000',time='2019-08-08 21:41')

-  **EphemJPL** – Class that controls the ephemeris querying in the JPL
   Horizons website.

   -  **\__init__()** – query the JPL Horizons website

      **INPUT**: *name*: the name of the object to search. *id_type*:
      the type of the name. For instance searching for ’Europa’, if the
      id_type is ’smallbody’ it will search for the asteroid 52 Europa,
      if the id_type is ’majorbody’ it will search for the Galilean
      satellite Europa. Default: ’majorbody’

   -  **get_position()** – Calculates and return the geocentric position
      of the target for the given time.

      **INPUT**: time (necessary).

      **OUTPUT**: Astropy SkyCoord Object with the geocentric coordinate
      of the object for the given times.

   -  **get_ksi_eta()** – returns the calculated :math:`\xi` and
      :math:`\eta` from input time, the star coordinate can be passed as
      well to skip the call of fit_d2_ksi_eta().

      **INPUT**: time (necessary), star coordinate (not necessary).

      **OUTPUT**: pair of :math:`\xi` and :math:`\eta`. If the input
      time is an array, two arrays for :math:`\xi` and :math:`\eta` will
      be returned.

   -  **get_pole_position_angle()** – Calculates and return the position
      and the aperture angle of the target’s pole for the given time.

      **INPUT**: pole: Pole coordinate (necessary); time: (necessary).

      **OUTPUT**: position angle of the pole in degrees; aperture angle
      in degrees.

.. code:: python

   from sora import EphemJPL
   from astropy.time import Time

   ephe3 = EphemJPL('Phoebe')

   ephe3.get_position('2019-06-07 03:54:22.60')

   ksi, eta = ephe3.get_ksi_eta(time='2019-06-07 03:54:22.60', star='19 21 18.63201 -21 44 25.3924')

   ephe3.get_pole_position_angle(pole='10 05 12.000 +41 28 48.000',time='2019-08-08 21:41')

.. _SubSubSec:code_classes_lightcurve:

Module: LightCurve
~~~~~~~~~~~~~~~~~~

-  **set_exposure** – Set the exposure time for the light curve as
   furnish by the user in seconds.

-  **set_filter** – Set the mean wavelength and width for the light
   curve as furnish by the user in microns.

-  **sigma_noise** – Set the :math:`\sigma_{LC}` as furnish by the user
   or calculat it based on the light curve flux.

-  **plot_lc** – Plots the observed light curve, and the best-fitted
   model if it was already determined.

-  **occ_detect** – Automatically search for the occultation event in
   the light curve using a BLS algorithm (Kovacs et al., 2002).

   **INPUTS**: *maximum_duration* (float, optional): Maximum duration of
   the occultation event (default is 1/4th of the light curve’s time
   span). *dur_step* (float, optional): Step size to sweep occultation
   duration event (default value is 1/2 of sampling). *snr_limit*
   (float,optional): Minimum occultation SNR. *n_detections* (int,
   optional): N best detections regardless from SNR. *n_detections* is
   superseded by *snr_limit*. **OUTPUT**: An ordered dictionary of
   :attr:‘name’::attr:‘value’ pairs for each Parameter.
   *occultation_duration*: time span of the occultation event.
   *central_time*: time of central instant of the event.
   *imersion_time*: instant of the begining of the event.
   *emersion_time*: instant of the end of the event. *time_err*:
   uncertainty of the event instants (from sampling). *depth*: depth
   magnitude of the occultation. *depth_err*: 1\ :math:`\sigma`
   uncertainty of the occultation depth. *baseline*: average flux level
   outside the occultation event. *baseline_err*: 1\ :math:`\sigma`
   uncertainty of baseline flux estimate. *snr*: Signal-to-noise ratio
   of the deep in respect with the baseline. *occ_mask*: Boolean mask
   with the lenght of the time series where data inside occultation is
   ‘True’. Note that when the multiple detection mode caution is
   required since other occultations in the data will influence
   statistic results.

   Single detection example:

   .. code:: python

      from sora import LightCurve
      import numpy as np

      input = 'LC_Chariklo_2017-04-09_Wabi.txt'
      time, arbtime, flux, dflux = np.loadtxt(input, usecols=(0,1,2,3), unpack=True)

      lc = LightCurve(time,flux)
      occ = lc.occ_detect()
      print(occ)
      #{'occultation_duration': 0.0004645648878067732,
      # 'central_time': 2457852.5916293273,
      # 'imersion_time': 2457852.5913970447,
      # 'emersion_time': 2457852.59186161,
      # 'time_err': 5.799811333417892e-07,
      # 'depth': 0.8663887801707082,
      # 'depth_err': 0.10972550419008305,
      # 'baseline': 0.9110181732552853,
      # 'baseline_err': 0.1904360360568157,
      # 'snr': 91.21719495827487,
      # 'occ_mask': array([False, False, False, ..., False, False, False])}

   Multiple detection using *snr_limit*:

   .. code:: python

      from sora import LightCurve
      import numpy as np

      input = 'LC_Chariklo_2017-04-09_Wabi.txt'
      time, arbtime, flux, dflux = np.loadtxt(input, usecols=(0,1,2,3), unpack=True)

      lc = LightCurve(time,flux)
      occ = lc.occ_detect(snr_limit=1)

      for i in range(len(occ['snr'])):
          print('\nDetection number: {}'.format(i+1))
          for key, value in occ.items():
              print('{}:{}'.format(key, value[i]))
      #Detection number: 1
      #occultation_duration:0.0004645648878067732
      #central_time:2457852.5916293273
      #imersion_time:2457852.5913970447
      #emersion_time:2457852.59186161
      #time_err:5.799811333417892e-07
      #depth:0.8663887801707082
      #depth_err:0.10972550419008305
      #baseline:0.9110181732552853
      #baseline_err:0.1904360360568157
      #snr:91.21719495827487
      #occ_mask:[False False False ... False False False]
      #
      #Detection number: 2
      #occultation_duration:1.1599622666835785e-05
      #central_time:2457852.592378953
      #imersion_time:2457852.592373153
      #emersion_time:2457852.592384753
      #time_err:5.799811333417892e-07
      #depth:0.472894269480926
      #depth_err:0.12751327397424012
      #baseline:0.8395567167536533
      #baseline_err:0.302469636172243
      #snr:5.185356378993047
      #occ_mask:[False False False ... False False False]
      #
      #Detection number: 3
      #occultation_duration:1.7399434000253677e-05
      #central_time:2457852.5908193258
      #imersion_time:2457852.590810626
      #emersion_time:2457852.5908280257
      #time_err:5.799811333417892e-07
      #depth:0.34605982894043885
      #depth_err:0.23829629773980976
      #baseline:0.8396265189404388
      #baseline_err:0.3025636864098951
      #snr:4.575034539625721
      #occ_mask:[False False False ... False False False]
      #
      #Detection number: 4
      #occultation_duration:0.0008438725490123034
      #central_time:2457852.591455855
      #imersion_time:2457852.591033919
      #emersion_time:2457852.591877791
      #time_err:5.799811333417892e-07
      #depth:0.49745658515149577
      #depth_err:0.43834959442481086
      #baseline:0.914112323959715
      #baseline_err:0.19028874519627717
      #snr:70.63232661080131
      #occ_mask:[False False False ... False False False]

-  **\_runBLS** – Private function used to find the best box fit
   suitable to the data.

-  **\_summarizeBLS** – Private function used to merge dictionaries
   returned by \_runBLS and keep values of common keys in list.

-  **occ_model** – Create a light curve occultation model considering
   the geometric box model, fresnel difraction, stellar diameter and
   exposure time.

   **INPUTS**: *self* (object): Object LightCurve. *t_ingress* (float):
   Ingrees time, in seconds. *t_egress* (float): Egress time, in
   seconds. *opa_ampli* (float) Opacity, opaque = 1.0, transparent =
   0.0. *mask* (array) Boolean mask with the lenght of the time series
   where data to be considered is ‘True’. *npt_star* (int) Number of
   subdivisions for computing the star size’s effects, default equal to
   12. *time_resolution_factor* (float) Steps for fresnel scale used for
   modelling the light curve, default equals to 10 steps for fresnel
   scale (or time exposure if smaller than fresnel scale). **OUTPUT**:
   *flux_inst* (array): Modelled Instrumental light flux. *time_model*
   (array): Modelled timing. *model_geometric* (array): Modelled light
   flux considering a box model. *flux_star* (array): Modelled light
   flux considering fresnel difraction and star’s diameter.
   *flux_fresnel* (array): Modelled light flux considering fresnel
   difraction.

   Example:

   .. code:: python

      from sora import LightCurve
      import numpy as np
      import matplotlib.pylab as pl

      inn = 'LC_Chariklo_2017-06-22_Hakos.txt'

      time, arbtime, flux, dflux = np.loadtxt(input, usecols=[0,1,2,3], unpack=True)

      lc = LightCurve(arbtime,flux)

      lc.exptime = 0.050                  # Exposure time in s
      lc.dist = 14.65925396               # Object distance in AU
      lc.vel  = 22.3572                   # Event velocity in km/s
      lc.d_star =  0.2104                 # Stelar radius in km

      t_ingress = 1270.4481
      t_egress = 1270.5985
      opacity = 0.435
      mask = (lc.time > t_ingress - 5) & (lc.time < t_egress + 5)

      lc.occ_model(t_ingress,t_egress,opacity)

      ##Plotting the observed and the modelled light curves, Figure 2a.
      ##Plotting the observed light curve, the box model, and the light curve affected by Fresnel diffraction, Figure 2b.

   |image| |image1|

-  **\__bar_fresnel** – Private function used to create the light curve
   affected by fresnel difraction.

-  **\__occ_model** – Private function used to create the light curve
   model during the fitting.

-  **occ_lcfit** – Function that fits the light curve parameters
   (ingress time, egress time and opacity) using a brute force
   minimisation of :math:`\chi^2`.

   **INPUTS**: *self* (object): Object LightCurve. *mask* (array)
   Boolean mask with the lenght of the time series where data to be
   considered is ‘True’. *t_ingress* (float): Initial guess for ingrees
   time, in seconds. *t_egress* (float): Initial guess for egress time,
   in seconds. *opa_ampli* (float) Initial guess for opacity, opaque =
   1.0, transparent = 0.0. *dt_ingress* (float): Bondary region for
   :math:`\chi^2` test, relative to the ingress time, in seconds.
   Default equals to zero means no variation *dt_egress* (float):
   Bondary region for :math:`\chi^2` test, relative to the egress time,
   in seconds. Default equals to zero means no variation *dopacity*
   (float) Bondary region for :math:`\chi^2` test, relative to the
   opacity. Default equals to zero means no variation

   **OUTPUT**: *chi2* (array): Tested chi squared values. *t_i* (array):
   Tested ingress values. *t_e* (array): Tested egress values. *opa*
   (array): Tested opacity values.

   Example only fitting Ingress time:

   .. code:: python

      from sora import LightCurve
      import numpy as np
      import matplotlib.pylab as pl

      inn = 'LC_Chariklo_2017-04-09_Weaver.txt'

      time, arbtime, flux, dflux = np.loadtxt(input, usecols=[0,1,2,3], unpack=True)

      lc = LightCurve(arbtime,flux)

      lc.exptime = 0.080                  # Exposure time in s
      lc.dist = 4.779391622236            # Object distance in AU
      lc.vel  = 15.4676223706             # Event velocity in km/s
      lc.d_star =  0.735                  # Stelar radius in km
      lc.sigma = 0.108

      #Only one side
      tmin = 7870
      tmax = 7890

      t_ingress = 7880.4792
      delta_t = 0.20

      chi2_ing =  lc.occ_lcfit(tmin=tmin, tmax=tmax, t_ingress=t_ingress, delta_t=delta_t)

      ##Plotting the Light Curve, Figure 3a
      ##Plotting the Chi squared curve for ingress time, Figure 3b

   |image2| |image3|

   Example of fitting all three parameters (Ingress time, egress time,
   and opacity):

   .. code:: python

      from sora import LightCurve
      import numpy as np
      import matplotlib.pylab as pl

      inn = 'LC_Chariklo_2017-04-09_Weaver.txt'

      time, arbtime, flux, dflux = np.loadtxt(input, usecols=[0,1,2,3], unpack=True)

      lc = LightCurve(arbtime,flux)

      lc.exptime = 0.080                  # Exposure time in s
      lc.dist = 4.779391622236            # Object distance in AU
      lc.vel  = 15.4676223706             # Event velocity in km/s
      lc.d_star =  0.735                  # Stelar radius in km
      lc.sigma = 0.108

      tmin = 7833
      tmax = 7845
      t_ingress = 7838.7
      t_egress = 7839.9
      opacity = 0.35

      delta_t  = 0.20 
      dopacity = 0.30 

      chi2_ring1 =  lc.occ_lcfit(tmin=tmin,tmax=tmax, t_ingress=t_ingress, t_egress=t_egress, opacity=opacity, delta_t=delta_t, dopacity=dopacity,loop=100000)

      ##Equivalent plots to the previous example
      ##Plotting the Light Curve, Figure 4a
      ##Plotting the Chi squared curve for ingress time, Figure 4b
      ##Plotting the Chi squared curve for egress time, Figure 4c
      ##Plotting the Chi squared curve for opacity, Figure 4d

   |image4| |image5| |image6| |image7|

.. _SubSubSec:code_classes_prediction:

Module: prediction
~~~~~~~~~~~~~~~~~~

-  Set of independent functions that are used by *Prediction* but do not
   depend on it.

   -  **occ_params()** – Calculates the parameters of the occultation.

      **INPUT**: *star*: The coordinate of the star in the same frame as
      the ephemeris. It must be a Star object. *ephem*: It must be an
      Ephemeris object.

      **OUTPUT**: instant of occultation, Closest Approach, Position
      Angle, geocentric relative velocity and object distance

      .. code:: python

         occ_params(star,eph)
         #(<Time object: scale='utc' format='jd' value=2458641.6600943254>,
         # <Quantity 0.16698067 arcsec>,
         # <Quantity 171.85352157 deg>,
         # <Quantity -16.85079849 km / s>,
         # <Distance 9.24102797 AU>)

   -  **prediction()**: Predicts occultations for the given inputs.

      **INPUT**: *ephem*: Ephemeris. It must be an EphemKernel object.
      *time_beg*: Initial time for prediction. *time_beg*: Final time
      for prediction. *mag_lim*: Faintest Gmag for search. *interv*:
      interval, in seconds, of ephem times for search. This is the step
      to generate ephemeris for fast search. *divs*: interval, in deg,
      for max search of stars. The ephemeris will be split in each “div”
      degrees to search occultations. *sigma*:Increase the range of
      distance for search based on sigma*ephem_error where ephem_error
      is the error given to EphemKernel object.

      **OUTPUT**: An Prediction object.

      .. code:: python

         from sora.prediction import prediction
         from sora.ephem import EphemKernel

         ephe2 = EphemKernel('Phoebe', 609, kernels=['path/ph19a.bsp', 'path/de438.bsp'])

         occs = prediction(ephe2, time_beg='2020-05-01 00:00:00.000', time_end='2020-10-01 00:00:00.000', mag_lim=17, divs=1)

         print(occs)
         #         Epoch             ICRS Star Coord at Epoch      ca  ...   G*     dist 
         #----------------------- ------------------------------ ----- ... ------ -------
         #2020-06-01 02:00:43.100 20 16 49.35703 -19 56 39.33936 0.054 ... 16.032   9.348
         #2020-06-11 17:02:49.950 20 15 02.29688 -20 03 38.58596 0.605 ... 15.324   9.231
         #2020-06-12 03:39:23.150 20 14 56.97040 -20 03 57.31622 0.670 ... 15.971   9.227
         #2020-07-24 12:57:53.650 20 02 57.42683 -20 44 41.96963 0.263 ... 13.371   9.042
         #2020-09-20 15:22:51.650 19 49 10.32921 -21 25 06.53792 0.154 ... 13.287   9.575

   -  **plot_occ_map()** – Plots an occultation map given input
      parameters.

      **INPUT**: *name*: Name of the object that is being plotted.
      *radius*: radius of the object, in km, to calculate size of path.
      *coord*: the GCRS coordinate of the star. *time*: The instant of
      closest approach, in the ISO format. *ca*: the closest approach
      distance, in arcsec. *pa*: the planet position angle with respect
      to the star at C/A, in deg. *vel*: the occultation velocity, in
      km/s. *dist*: the planet distance, in AU. *\**kwargs*: other
      parameters for plot control.

      Please refer to to a full example of all kwargs present in
      *plot_occ_map()*

-  **PredictionTable**

   -  **\__init__()** – Instantiate PredictionTable object.

      **INPUT**: *time*: the instant of geocentric closest approach;
      *coord_star*: the coordinate of the star; *coord_obj*: the
      geocentric coordinate of the object, *ca*: geocentric closest
      approach distance; *pa*: position angle at ca; *vel*: occultation
      geocentric velocity; *dist*: Object distance, in AU; *mag* and
      *mag_20* for the magnitude and magnitude normalised to 20km/s. At
      least one of these parameters are required; *long*: Longitude of
      the sub-star point at ca (optional); *loct*: Local solar time at
      long (optional). *\**kwargs*: any other parameter will be parsed
      directly by Astropy Table.

      **The Prediction Table is not intended to be directly instantiated
      by the user.**

   -  **to_praia()** – Saves the table in the PRAIA format

      **INPUT**: *filename*: Name of the file to save table

   -  **from_praia()** – Creates a Prediction table from PRAIA
      prediction table.

      **INPUT**: *filename*: Name of the file in PRAIA format; *name*:
      The name of the object predicted; *radius*: the radius of the
      object (optional).

      **OUTPUT**: A Prediction table.

      .. code:: python

         from sora.prediction import PredictionTable

         occs = PredictionTable.from_praia('g4_occ_datauc4_HIM_G2_table', name='Himalia', radius=85)
         print(occs)
         #         Epoch             ICRS Star Coord at Epoch    ... long  loct
         #----------------------- ------------------------------ ... ---- -----
         #2020-04-06 02:05:35.000 19 46 15.30170 -21 40 23.57500 ...   70 06:46
         #2020-04-08 06:19:54.000 19 47 25.29390 -21 38 18.27800 ...    5 06:39
         #2020-04-14 07:46:30.000 19 50 24.76210 -21 32 46.42700 ...  338 06:18
         #2020-04-27 13:53:34.000 19 55 30.85740 -21 22 31.22900 ...  234 05:31

   -  **to_ow()** – Creates the files ‘tableOccult_update.txt‘ and
      ‘LOG.dat‘ to feed Occult Watcher.

      **INPUT**: *ow_des*: the OW designation of the object; *mode*:
      “append” if the output must append existing files, “restart” to
      overwrite existing files. Default: “append”.

   -  **plot_occ_map()** – Plots the maps for all the occultations. All
      the parameters required by *plot_occ_map()* are automatically
      filled. The user can pass the parameters for plot control. Please
      refer to Section `5.5.8 <#SubSubSec:code_maps>`__ to a full
      example of all kwargs present in *plot_occ_map()*

      .. code:: python

         from sora.prediction import PredictionTable

         occs = PredictionTable.from_praia('g4_occ_datauc4_HIM_G2_table', name='Himalia', radius=85)
         occs.plot_occ_map()
         ## This will plot all the occultations in the occs.

         occs[0].plot_occ_map()
         ## This will plot only the first occultation in occ

         occs['2020-04-06'].plot_occ_map()
         ## This will plot all the occultations predicted for April 06, 2020.

         occs['2020-04'].plot_occ_map()
         ## This will plot all the occultations predicted for April 2020.

         occs['2020-04-06 02'].plot_occ_map()
         ## This will plot all the occultations predicted for April 06, 2020 between 02 and 03 UT.

   -  **other methods** – Prediction inherits all attributes and methods
      from Astropy.table.Table. It is expected that all methods from
      Table should work fine. This includes the function to write table
      in other formats, such as *csv* and *latex* with **write()**,
      pretty printing with **pprint_all()**, stacking Prediction tables
      with **vstack()**, etc. Please refer to
      https://docs.astropy.org/en/stable/table/ for all Table
      Documentation.

.. _SubSubSec:code_classes_occultation:

Module: occ_fit
~~~~~~~~~~~~~~~

-  Set of independent functions that are used by *Occultation* but do
   not depend on it.

   -  **positionv()** – Calculates the position and velocity of the
      occultation shadow relative to the observer.

      **INPUT**: *star*: The coordinate of the star in the same frame as
      the ephemeris. It must be a Star object. *ephem*: It must be an
      Ephemeris object. *observer* (Observer): It must be an Observer
      object. *time*: Instant to calculate position and velocity

      **OUTPUT**: The orthographic projection of the shadow relative to
      the observer (f and g) and velocity in these directions.

      .. code:: python

         from sora import Observer, Star, Ephemeris, positionv

         casleo = Observer('Casleo', '-69 17 44.9', '-31 47 55.6', 2492)
         star = Star(coord='19 21 18.63201 -21 44 25.3924')
         eph = Ephemeris('Phoebe', 'ephem_phoebe_ph15_jd.txt')

         f,g,vf,vg = positionv(star=star,ephem=eph, observer=casleo, time='2019-06-07 03:54:22.60')
         print(f,g,vf,vg)
         #-225.359 -48.573 16.957 2.487

   -  **fit_ellipse()** – fit ellipse to given points.

      **INPUT**: list of Occultation objects. The position will be read
      automaticaly from these objects.

      Required parameters: *center_f*, *center_g*, *equatorial_radius*,
      *oblateness*, *pos_angle* as initial parameters. For each of
      theses parameters, the interval of values must be given by
      *dcenter_f*, *dcenter_g*, *dequatorial_radius*, *doblateness*,
      *dpos_angle*. If any delta is not given, it will be set to zero.

      *loop*: Number of different ellipses to try. Default: 10,000,000

      *dchi_min*: if a number is set to this parameter, fit_ellipse will
      only saves values where
      :math:`\chi^2 < \chi^2 + \textrm{dchi\_min}`. In this case, it
      will repeatedly generate random ellipses until *number_chi* is
      reached. *number_chi*: minimum number of ellipsis within
      *dchi_min*. Default: 10,000.

      *log*: if True, it will print the steps of the fitting. Default:
      False.

      **OUTPUT**: a *ChiSquare* object.

      .. code:: python

         from sora.occ_fit import fit_ellipse, Occultation

         occ1 = Occultation(star1, ephem1)
         # add observations for occ1
         occ2 = Occultation(star2, ephem2)
         # add observationn for occ2

         chi2  = fit_ellipse(occ1, occ2, center_f=87, center_g=-60, dcenter_f=3, dcenter_g=3, equatorial_radius=84, dequatorial_radius=3,oblateness=0.23, doblateness=0.1, pos_angle=-5, dpos_angle=10 ,loop=10000000, dchi_min=10)

   -  **\_PositionDict(dict)** – A modified dict class to handle
      positions. It should not be used by the user.

-  **Occultation**

   -  **\__init__()** – instantiate the star object.

      **INPUT**: *star*: The coordinate of the star in the same frame as
      the ephemeris. It must be a Star object. *ephem*: Ephemeris. It
      must be one of the Ephemeris objects.

   -  **add_observation()** – Add Observers to the Occultation object.

      **INPUT**: *obs*: It must be an Observer object; *LightCurve*: It
      must be an LightCurve object

   -  **observations()** – Shows all observations added to the
      occultation. No input required.

   -  **remove_observation()** – Remove an observation from the
      occultation list.

      **INPUT**: *key*: The name given to Observer or LightCurve to
      remove from the list; *key_lc*: In the case where repeated names
      are present for different observations, *key_lc* must be given for
      the name of the LightCurve and *key* will be used for the name of
      the Observer.

   -  **fit_ellipse()** – Calls the *fit_ellipse* function using only
      this occultation. Any input parameter will be parsed directly by
      *fit_ellipse*.

   -  **new_astrometric_position()** – Calculates and print new
      astrometric position given the params fitted by *fit_ellipse*.

      **INPUT**: *time*: the time the user wants to the determine the
      position. Default: Time at Closest Approach; *offset*: The offset
      the user wants to give to the ephemeris. If not given, it will use
      the offset calculated by *fit_ellipse*. If no fit have been made,
      it uses zero. The offset param must be a list with 3 values and
      can be in distance units for the offset in the projection, example
      “offset=[f, g, ’km’]” where X and Y are in km in the direction f
      (norte) and g (east). It can also be an angle value for angular
      offsets, example “offset=[da_cos_dec, d_dec, ’mas’]”; *error*: the
      error for the offset applied. Its use is similar to the *offset*
      param.

   -  **check_velocities()** – Shows the radial velocities for the
      immersion and emersion for each positive chord given the center
      calculated by fit_ellipse.

      **INPUT**: None.

   -  **plot_chords()** – Makes the plot of the chords of the
      observations added to the Occultation.

      **INPUT**: *all_chords*: Default=True, if it is to print all the
      chords. If false, it does not plot the points or chords disabled
      by the user; *positive_color*: The color to plot positive colors.
      Default: ’blue’; *negative_color*: The color to plot negative
      colors. Default: ’green’; *error_color*: The color to plot the
      error bars of the points. Default: ’red’.

   -  **positions** – This is a function that works as dictionary. Once
      called, it will calculated the position and velocity for all
      immersions and emersions for the light curves added to the
      Occultation. It uses the *\_PositionDict* object. It is used to
      able or disable points for the fitting process.

.. code:: python

   from sora import Observer, Star, EphemPlanete, Occultation, LightCurve

   star = Star(coord='19 21 18.63201 -21 44 25.3924')
   eph = Ephemeris('Phoebe', 'ephem_phoebe_ph15_jd.txt')

   casleo = Observer('Casleo', '-69 17 44.9', '-31 47 55.6', 2492)
   casleo_lc = LightCurve('Casleo', time, flux_input)

   occ = Occultation(star=star, ephem=eph)
   occ.add_observation(casleo, casleo_lc)

   # to disable casleo observation to not use it in the fitting process.
   occ.positions['Casleo'] = 'off'
   # instead, to disable only the immersion of casleo observation.
   occ.positions['Casleo']['immersion'] = 'off'

.. _SubSubSec:code_maps:

Map documentation
~~~~~~~~~~~~~~~~~

To generate prediction or post-fit occultation maps, the function to be
used is *sora.prediction.plot_occ_map(*). This function receives many
required parameters that states the orientation and path of the
occultation. When using this function within sora.occ_fit.plot_occ_map()
or sora.Prediction.plot_occ_map() all the required parameters are
automatically send to the original function. But the user is still able
to pass many kwargs are used to configure the occultation map. Without
any configuration parameter, the map will have the plain view of the
Earth and path. For example.

.. code:: python

   ### from a Prediction Table of Phoebe
   occs.plot_occ_map()
   ## Phoebe_2018-06-19T04:36:56.400.png generated

|image8|

-  nameimg: Change the name of the imaged saved.

   .. code:: python

      occs.plot_occ_map(nameimg='Phoebe_teste')
      ## Phoebe_teste.png generated

-  resolution: Identify the type of cartopy feature resolution. "1"
   means a resolution of "10m", "2" a resolution of "50m" and "3" a
   resolution of "100m". The default is "2". .

   .. code:: python

      occs.plot_occ_map(resolution=1)
      occs.plot_occ_map(resolution=3)

   |image9| |image10|

-  states: Plots the state division of the countries. The states of some
   countries will only be shown depending on the resolution. For
   instance, USA states are shown for all resolutions, but Brazilian
   states will be shown only for resolution equal to "1" or "2". This is
   a cartopy characteristics and cannot be changed. Default=True.

   .. code:: python

      occs.plot_occ_map(states=True)
      occs.plot_occ_map(states=False)

-  zoom: Zooms in or out of the map. It must be a number. For the number
   given in zoom, the dimensions of the map will be divided by that
   number. For instance, if zoom=2, it will be shown the Earth divided
   by 2 in X and Y. Default=1. .

   .. code:: python

      occs.plot_occ_map(zoom=2)
      occs.plot_occ_map(zoom=0.5)

   |image11| |image12|

-  centermap_geo: Center the map given coordinates in longitude and
   latitude. It must be a list with two numbers. Default=None. If
   coordinate is on the other side of the map, it gives an error.
   (left).

-  centermap_delta: Displace the center of the map given displacement in
   X and Y, in km. It must be a list with two numbers. Default=None.
   (right).

   centermap_geo and centermap_delta are only a displacement on original
   projection. If Earth rotation is needed, please see centerproj.

   .. code:: python

      occs.plot_occ_map(centermap_geo=[0,-40], zoom=2)
      # Center the map on longitude=0.0 deg and latitude=-40 deg.
      occs.plot_occ_map(centermap_delta=[5000,-1000], zoom=2)
      # Displace the center of the map 5000 km East and 1000 km South

   |image13| |image14|

-  centerproj: Rotates the Earth to show occultation with the center
   projected at a given longitude and latitude. (left).

-  labels: Plots text above and below the map with the occultation
   parameters. Default=True. (right).

-  meridians and parallels: Plots lines representing the meridians and
   parallels given such interval. Default=30 for both parameters. So it
   will plot lines representing these values each 30º. (right).

   .. code:: python

      occs.plot_occ_map(centerproj=[0,-40])
      # Rotate to center the map projection on longitude=0.0 deg and latitude=-40 deg.
      occs.plot_occ_map(labels=False, meridian=10, parallels=10, zoom=2)
      # Displace the center of the map 5000 km East and 1000 km South

   |image15| |image16|

-  sites: Plots site positions in map. It must be a python dictionary
   where the key is the name of the site, and the value is a list with
   longitude, latitude, delta_x, delta_y and color. *delta_x* and
   *delta_y* are displacement, in km, from the point of the site in the
   map and the name. *color* is the color of the point. (left).

-  countries: Plots the names of countries. It must be a python
   dictionary where the key is the name of the country and the value is
   a list with longitude and latitude of the lower left part of the
   text. (right).

   .. code:: python

      sites = {}
      sites['Foz'] = [ -54.5936, -25.4347, 10, 10, 'blue']
      sites['SOAR'] = [ -70.73919, -30.238027, 10,10,'green']
      sites['La Silla'] = [-70.7393888, -29.254611, 10,10,'blue']
      occs.plot_occ_map(zoom=5, labels=False, sites=sites)

      countries = {}
      countries['Brazil'] = [-52.5983973, -23.5570511]
      countries['Argentina'] = [-67.2088692, -35.1237852]
      occs.plot_occ_map(zoom=3, labels=False, countries=countries, states=False)

   |image17| |image18|

-  offset: applies an offset to the ephemeris, calculating new CA and
   instant of CA. It is a pair of delta_RA*cosDEC and delta_DEC. (left).

-  mapstyle: Define the color style of the map. 1 is the default black
   and white scale. 2 is a colored map. (right).

   .. code:: python

      occs.plot_occ_map(zoom=3, offset=[-40,50])
      # Applies offsets of -40 mas in Delta_alpha_cos_delta and 50 mas in Delta_delta
      occs.plot_occ_map(zoom=3, mapstyle=2)
      # Plots a colored map, without offset

   |image19| |image20|

-  error: Ephemeris error in mas. It plots a dashed line representing
   radius + error. To change the color of these lines, the name of the
   color must be given to lncolor. (left).

-  ring: Similarly to error, it plots a dashed line representing the
   location of a ring. It is given in km, from the center. To change the
   color of these lines, the name of the color must be given to rncolor.

-  atm: Similarly to error, it plots a dashed line representing the
   limit of an atmosphere. It is given in km, from the center. To change
   the color of these lines, the name of the color must be given to
   atmcolor.

-  heights: It plots a circular dashed line showing the locations where
   the observer would observe the occultation at a given height above
   the horizons. This must be a list. To change the color of these
   lines, the name of the color must be given to hcolor. (right).

   .. code:: python

      occs.plot_occ_map(zoom=3, labels=False, error=15)
      # Shows an error bar of 15 mas
      occs.plot_occ_map(heights=[30])
      # Shows where the observer will see the occultation with a 30deg height above the horizons.

   |image21| |image22|

-  mapsize: The size of figure, in cm. It must be a list with two
   values. Default = [46.0, 38.0].

-  cpoints: Interval for the small points marking the center of shadow,
   in seconds. Default=60. To change the color of these points, the name
   of the color must be given to ptcolor.

-  alpha: The transparency of the night shade, where 0.0 is full
   transparency and 1.0 is full black. Default = 0.2.

-  fmt: The format to save the image. It is parsed directly by
   matplotlib.pyplot. Default = ’png’.

-  dpi: "Dots per inch". It defines the quality of the image. Default =
   100.

-  nscale, cscale, sscale and pscale: Arbitrary scale for the size for
   the name of the site, for the name of the country, for the size of
   point of the site, and the size of the points that represent the
   center of the shadow, respectively. This scale is arbitrary and is
   proportional to the size of the image.

-  lncolor, outcolor: To change the color of the line that represents
   the limits of the shadow over Earth and the color of the lines that
   represents the limits of the shadow outside Earth, respectively.

Python hints, tips and tools (for users) [Sec:Python_hints]
===========================================================

Here we will describe some hints ant tips for python usage with SORA.
Please feel free to send us any T&H you thing it is usefull.

Reporting a bur or issue [SubSec:tips_report_bug]
-------------------------------------------------

The user has a few different ways to report problems and bugs in the
code. Preferably using GITm but one can either report via email to the
developers or via SLACK (for LIneA members) on channel #sora. Using GIT
follow the steps:

#. Login into your GIT account on https://github.com/

#. Change to the SORA repository

#. Click on “Issues” and then on the button “New Issue”.

#. Put a title and a description of a problem, bug or error you are
   facing with SORA and submit.

#. The developers team will work on a solution and the fix will be
   available in a future release.

Jupyter Notebooks [SubSec:tips_jupyter_notebooks]
-------------------------------------------------

Jupyter Notebook is a web-based interactive computational environment.
They are made of an ordered list of input/output cells which can contain
code, text (using Markdown), plots etc.

Jupyter Notebook is installed by default if you use Anaconda. However,
JupyterLab may not.

.. code:: python

   $ conda install -c conda-forge jupyterlab

Alternatively, using pip:

.. code:: python

   $ pip install jupyterlab

Once you have it installed it can be launched as follows

.. code:: python

   # Jupyter Notebook:
   $ jupyter notebook
   # or lauching a specific file locally:
   $ jupyter notebook file.ipynb

   # JupyterLab:
   $ jupyter lab

JupyterLab provides flexible building blocks for interactive,
exploratory computing, compared to the jupyter notebook interface:
https://jupyterlab.readthedocs.io/en/stable/user/interface.html

Some shortcuts:

-  Run cell: shift + enter

-  Run cell and create a new below: ctrl + enter

-  Change to Markdown: m

-  Add cell above: a

-  Delete cell: d + d

-  Add cell below: b

.. [1]
   References given here are the ones up to the first version of SORA on
   May/2020. The international collaboration published many other works
   after that.

.. [2]
   Born & Wolf 1980 -
   https://ui.adsabs.harvard.edu/abs/1980poet.book.....B/abstract

.. [3]
   Gaia Collaboration 2016, 2016a and 2018. Mission:
   https://ui.adsabs.harvard.edu/abs/2016A%26A...595A...1G/abstract;
   DR1:
   https://ui.adsabs.harvard.edu/abs/2016A%26A...595A...2G/abstract;
   DR2: https://ui.adsabs.harvard.edu/abs/2018A%26A...616A...1G/abstract

.. [4]
   Zacharias et al. 2004 -
   https://ui.adsabs.harvard.edu/abs/2004AAS...205.4815Z/abstract

.. [5]
   van Belle 1999 -
   https://ui.adsabs.harvard.edu/abs/1999PASP..111.1515V/abstract

.. [6]
   Kervella et al. 2004 -
   https://ui.adsabs.harvard.edu/abs/2004A%26A...426..297K/abstract

.. [7]
   In Windows systems, once the user installs Anaconda, it is possible
   to access a unix-like terminal and do this procedure.

.. [8]
   When facing any issue, the user can try installing each package
   individually, instead of all packages in the same command line.

.. |image| image:: static/img/exemplo_lcmodel_1.png
.. |image1| image:: static/img/exemplo_lcmodel_2.png
.. |image2| image:: static/img/exemplo_lcfit_1a.png
.. |image3| image:: static/img/exemplo_lcfit_1b.png
.. |image4| image:: static/img/exemplo_lcfit_2a.png
.. |image5| image:: static/img/exemplo_lcfit_2b.png
.. |image6| image:: static/img/exemplo_lcfit_2c.png
.. |image7| image:: static/img/exemplo_lcfit_2d.png
.. |image8| image:: static/img/maps/Phoebe_plain.png
.. |image9| image:: static/img/maps/Phoebe_res1.png
   :width: 49.0%
.. |image10| image:: static/img/maps/Phoebe_res3.png
   :width: 49.0%
.. |image11| image:: static/img/maps/Phoebe_zoom2.png
   :width: 49.0%
.. |image12| image:: static/img/maps/Phoebe_zoom05.png
   :width: 49.0%
.. |image13| image:: static/img/maps/Phoebe_center_geo.png
   :width: 49.0%
.. |image14| image:: static/img/maps/Phoebe_center_delta.png
   :width: 49.0%
.. |image15| image:: static/img/maps/Phoebe_centerproj.png
   :width: 49.0%
.. |image16| image:: static/img/maps/Phoebe_merpar.png
   :width: 49.0%
.. |image17| image:: static/img/maps/Phoebe_sites.png
   :width: 49.0%2
.. |image18| image:: static/img/maps/Phoebe_country.png
   :width: 49.0%
.. |image19| image:: static/img/maps/Phoebe_offset.png
   :width: 49.0%
.. |image20| image:: static/img/maps/Phoebe_mapstyle.png
   :width: 49.0%
.. |image21| image:: static/img/maps/Phoebe_error.png
   :width: 49.0%
.. |image22| image:: static/img/maps/Phoebe_height.png
   :width: 49.0%

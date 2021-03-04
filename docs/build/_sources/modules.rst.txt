.. _Sec:modules:

Modules
=====================

To perform the complete data reduction, pre- and post-occultation, we
have implemented the modules described below. Each of them have their
own set of input and output values (I/O). It is important to note that
those I/O can be the same for different modules or even some outputs
from a module will work as input for other, meaning that they are
integrated and has an inside communication. The modules are divided as
follows:

-  **Implemented modules:**

   #. **Prediction:** Module to make the prediction of stellar
      occultation. It uses the object ephemeris (plus its updates and
      offsets) and a star catalogue to create prediction maps of the
      candidate events. The maps present the object shadow path on Earth
      and contain the occultation related information.

   #. **Light Curve (LC fit):** Module to determine the ingress and
      egress instants from an observer data-set (the light curve). It
      uses the light curve (with time and normalised relative flux
      between the target star and one or more reference stars) to
      calculate the instants that the object enters in front of the star
      and leave (immersion and emersion instants, respectively).

   #. **Sky Projection:** Module to project each of the times obtained
      with the module *LC fit* into positions (:math:`f_i` and
      :math:`g_i`) in the sky plane for a Geocentric reference frame. It
      uses the the times obtained in module *LC fit* and the occultation
      information from module *Prediction*, added to observer
      information, to convert the instants into positions in the sky
      plane

   #. **Ellipse fit:** Module to fit an ellipse with given points. It
      uses the projected positions in the sky plane (:math:`f_i` and
      :math:`g_i`) obtained in module *Sky Projection* to fit the 5
      parameters of one ellipse. (See specific section to know more
      about it, e.g., what assumptions are made if less than 5 points
      are available).

In the next chapters we present a description of each of the modules and
the algorithm and processes done in each part. Note that SORA is mostly
based in the processes and science from Bruno Sicardy’s programs, but it
is not a “literal translation” from fortran to python, i.e., SORA has
its own routines, features, and other approaches in the methods to
calculate the outputs.

Prediction [SubSec:prediction]
------------------------------

Module to make the prediction of stellar occultations. It uses the
object ephemerid (plus its updates and offsets) and a star catalogue
(Gaia-DR2) to calculate the occultation parameters and create the
prediction maps of the candidate events.

The ephemerid is calculated for a time interval provided by the user
using the ephemerid kernel (bsp files), or through a query to
Horizons/JPL. Those positions will be projected in the sky plane for
each instant given. The projected positions can also be calculated using
the equations given in Bruno’s Sicardy’s fortran code ephem_planete,
which uses a seconde degree polynom within 2 hours centered in the event
time.

The star information is obtained from Gaia DR2 and NOMAD catalogues
through a query to
`VizieR <http://vizier.u-strasbg.fr/viz-bin/VizieR>`__ database. The
outputs are the occultation parameters and the prediction maps, which
present the object shadow path on Earth and contain all the occultation
parameters (besides object information some of the parameters are:
occultation date, star RA and DEC, C/A, P/A, event relative velocity,
object distance in au, normalised magnitude).

This module can also import files and tables in PRAIA format and can
create output files in the user preference:

-  tables in PRAIA format;

-  files to feed Occult Watcher;

Kernels and ephemerids [SubSec:Kernels_bsp]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SORA can use the ephemerid file the user obtained from a query to
Horizons/JPL or can use the ephemerid kernels (bsp files). In this last
case, the user needs to have the files locally (i.e., in the same
environment: personal computer, online jupyter-notebook tool, etc.).

The BSP files from JPL can be found at:

-  Planetary ephemedides:

   https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/

-  Satellite ephemerides:

   https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/

It is important to read the “\*.cmt” files to know which objects are
contained within the BSP file. The irregular satellites of Jupiter, for
example, are within JUP343, but the inner satellites are within JUP310.
Also, note that a higher number does not mean that it is the most up to
date to all satellites.

The BSP files from NIMA - Lucky Star (by Josselin Desmars) can be found
at:

-  Planetary ephemerides (mostly minor bodies):

   https://lesia.obspm.fr/lucky-star/astrom.php

   https://lesia.obspm.fr/lucky-star/nima.php

-  Phoebe (Saturn IX):

   http://josselin.desmars.free.fr/work/research/ph20/

In those address, the user just need to select the object and in the
individual page, go to “Ephemeris & Observations” section and find the
“bsp file”.

The BSP files from INPOP (IMCCE) can be found at:

-  Planetary ephemerides:

   https://www.imcce.fr/recherche/equipes/asd/inpop/download19a

-  Satellite ephemerides:

   http://nsdb.imcce.fr/multisat/nssephmf.htm

Note that in the case of INPOP there is no BSP. The user will have to
generate the ephemeride file (using other tools) and then use it on SORA
(it is important to notice that “EphemPlanete” does not work for
prediction, only for occultation)

Other ephemerid files (such as NOE from Valéry Lenay) should be asked
directly to the author.

SPKID (object code) [SubSec:spkid]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to identify the object code (SPKID) for the minor bodies, the
user can search on the following address:
https://ssd.jpl.nasa.gov/sbdb.cgi .

For planets and Pluto the SPKID is the order from the Sun \*100 + 99.
Examples:

-  Mercury = 199

-  Venus = 299

-  Earth = 399

For the planet’s satellites the SPKID is the order of the planet from
the Sun \*100 + the number of the object.

-  Himalia (Jupiter VI) = 506

-  Phoebe (Saturn IX) = 609

-  Charon (Pluto I ) = 901

Light Curve fit - LC fit [SubSec:LC_fit]
----------------------------------------

This module intends to determine the ingress and egress instants from an
observer data-set (the light curve). It uses the time (usually seconds
with respect to 00:00 UTC from the occultation day), and the normalised
relative flux between the target star (plus body) and one or more
reference stars, to calculate a square well model convoluted with the
effect of star diam (diffraction), Fresnel scale, and integration time.
If the flux given is not normalised or it has some trend (due to sky
variations, for example), the user can choose to normalise the light
curve with a polynomial function (using the function normalize). The
synthetic light curve is then compared with the observation to obtain
the instants that the object enters in front of the star and leave
(immersion and emersion instants, respectively). This is an automated
and updated software based on Bruno Sicardy’s program bar.f (and its
variations).

Before considering any effect, the program uses a simple square well
function to search for significant drops in the light curve. The user
can specify the number of drops the program will look for and then the
program will consider each drop separately to convolute with the effects
that can affect the light curve. In the version 0.1, only the most
significant drop in the light curve will be considered to create a
synthetic light curve.

The first effect to be considered in the synthetic light curve is the
Fresnel diffraction scale (or simply Fresnel scale). The light emitted
by a point source (assumed to be at the infinity to yield planar waves)
incident on a sharp-edged obstacle (such as a TNO) is diffracted. Owing
to the Huygens-Fresnel principle of wave propagation, each point of a
wave front may be considered as the center of a secondary disturbance
giving rise to spherical wavelets, which mutually interfere. If part of
the original wave front is blocked by an obstacle, the system of
secondary waves is incomplete, so that diffraction occurs. When observed
at a finite distance D from the obstacle, this effect is known as
“Fresnel diffraction”, and falls within the scope of the Kirchhoff
diffraction theory which remains valid as long as the dimensions of the
diffracting obstacles are large compared to the observed wavelength
:math:`\lambda` and small compared to D [2]_. The characteristic scale
of the Fresnel diffraction effect is
:math:`\sqrt{\frac{\lambda * D}{2}}` (where :math:`\lambda` is the
wavelength of the observation, and :math:`D` is the geocentric distance
of the object). The Fresnel scale at 40 AU, observed at
:math:`\lambda = 0.4 \mu`\ m, is :math:`\sim`\ 1.1 km, so the
diffraction must be seriously taken into account to analyse the
occultations by TNOs.

In this routine, SORA calculates automatically the Fresnel scale for
wavelengths of 0.4 :math:`microns` and 1.0 :math:`microns`
(corresponding to 0.7 :math:`\pm` 0.3 :math:`microns`) – if the user do
not provide other values – to create two synthetic light curves, and
then get the average value for each point in the LC.

The second effect to be considered in the light curve is the star
angular diameter, which is the star size projected at the object
distance, as seen by an observer on Earth. This value is calculated on
nodule **Sky Projection**
(sec. `[SubSec:sky_projection] <#SubSec:sky_projection>`__). In this
routine, for each observational point, we divide the star size (assumed
as an uniform sphere) in 24 fractions and check how much of it is
occulted by the object, here assumed to be a sphere.

The last effect involves the exposure time together with the event
velocity. The synthetic curve is created with a time resolution that is
proportional to 1/10 of the Fresnel scale divided by the event velocity,
calculated in the prediction. The final effect for each observational
point is then considering the average of the number of synthetic points
within the exposure time (not considering cycle or dead times).

The final synthetic light curve is the convolution of those three
effects and it is compared with the observational light curve using a
:math:`\chi^2` minimisation method. The minimum :math:`\chi^2` gives the
best fit for the immersion and emersion times, with their uncertainties.

Sky Projection [SubSec:sky_projection]
--------------------------------------

This module is an automated and updated software based on Bruno
Sicardy’s programs (positionv.f, ephem_planete.f, fit_d2_ksi_eta.f,
polyfit.f) and Julio Camargo’s gstar.f. It calculates the projected
position in the sky plane of the object limb relative to the star
position, considering the times of immersion and emersion obtained in
*LC fit*, and using the geocenter as reference frame.

This module is divided in many different small routines (created as
python objects: `[item:routine:ephem] <#item:routine:ephem>`__ – Ephem;
`[item:routine:star_position] <#item:routine:star_position>`__ – Star
position;
`[item:routine:sky_projection] <#item:routine:sky_projection>`__ – Sky
Projection; `[item:routine:star_diam] <#item:routine:star_diam>`__ –
Star Diam;

#. **Ephem** [item:routine:ephem]

   This routine generates geocentric positions for the object around the
   event time. It also generates geocentric position for the Sun at the
   event predicted time.

#. **Star position** [item:routine:star_position]

   This routine calculates the star position on the occultation date,
   propagated from the catalogue epoch (e.g., epoch for GDR2 is 2015.5)
   to the event epoch, considering star parameters such as proper motion
   and parallax.

#. **Sky Projection** [item:routine:sky_projection]

   This routine uses the object ephemerid (all object positions -
   RA,DEC), and project in the sky plane to obtain the positions
   (:math:`\xi, \eta`) using the equations below:

   .. math::

      \centering
          \label{Eq:xi_projected}
          \xi = -D*cos(\delta_{obj})*sin(\alpha_{obj}-\alpha_{star})

   .. math::

      \centering
          \label{Eq:eta_projected}
          \eta = -D*[sin(\delta_{obj}-\delta_{star}) + 2*cos(\delta_{obj})*sin(\delta_{star})*sin^2(
          \frac{\alpha_{obj}-\alpha_{star}}{2})

   Then the routine will fit a second degree polonium to the apparent
   path in the sky plane (:math:`a \cdot x^2 + b \cdot x + c`) and give
   the coefficients :math:`a_\xi , b_\xi , c_\xi` and
   :math:`a_\eta , b_\eta ,c_\eta`.

#. **Star Diam** [item:routine:star_diam]

   This routine is intended to calculate diameter of the star projected
   at object distance using three methods. Preferentially the star
   diameter is calculated using Gaia-DR2 [3]_ information (when
   available). Also it uses the star magnitudes B, V, and K obtained
   from NOMAD catalogue [4]_ and equations from Van Belle (1999) [5]_
   and Kervella et. al (2004) [6]_.

   #. Van Belle (1999) [item:star_diam_vanBelle]

      This method uses equations from van Belle (1999) – Publi. Astron.
      Soc. Pacific 111, 1515-1523:

      .. math::

         \centering
             \label{Eq:van_belle}
         \begin{split}
             Diam_V = 10^{A_V + B_V*(V - K) -0.2*V} \\
             Diam_B = 10^{A_B + B_B*(B - K) -0.2*B}
         \end{split}

      where B, V and K are the star magnitudes in those bands and

      -  For Super Giant star: :math:`A_V`\ = 0.669 , :math:`B_V`\ =
         0.223, :math:`A_B`\ = 0.648, :math:`B_B`\ = 0.220

      -  For Main Sequence star::math:`A_V`\ = 0.500 , :math:`B_V`\ =
         0.264, :math:`A_B`\ = 0.500, :math:`B_B`\ = 0.290

      -  For Variable star::math:`A_V`\ = 0.789 , :math:`B_V`\ = 0.218,
         :math:`A_B`\ = 0.840, :math:`B_B`\ = 0.211

      The routine calculate the star diameter assuming a super giant to
      overestimate the star size. This is because the size of the star
      will be usually negligible compared to other effects (Fresnel
      scale or exposure time). If the star size is big enough to get
      some influence on the light curve, observation of the star and a
      better determination of its size will be needed.

   #. Kervella (2004) [item:star_diam_Kervella]

      This method uses equations from Kervella et al. (2004) – A&A Vol.
      426, No. 1 to calculate the star diameter:

      .. math::

         \centering
             \label{Eq:kervella}
         \begin{split}
             Diam_V = 10^{0.0755*(V − K) + 0.5170 − 0.2*K} \\
             Diam_B = 10^{0.0535*(B − K) + 0.5159 − 0.2*K} 
         \end{split}

      where B, V and K are the star magnitudes in those bands.

   #. Gaia (2018) [item:star_diam_Gaia]

      This method uses information from the GDR2 catalogue. It will give
      the star radius and distance (from parallax) and using the
      effective temperature :math:`T_{eff}` calculated using
      *Apsis-Priam*, the Radius and Luminosity calculated using
      *Apsis-Flame* (Refs: Andrae et al. A&A 616, A8 (2018);
      Bailer-Jones et al. A&A 559 (2018); Bailer-Jones et al. AJ, V.
      156, Issue 2, id. 58, 11pp (2018)) it gives de star diameter as:

      .. math::

         \centering
             \label{Eq:gaia}
             Diam = 2 * arctan \left(\frac{R_*}{2 D}\right)

      where R\ :math:`_*` is the star radius in solar radius and D the
      object distance in km.

#. **Positionv** [item:routine:positionv]

   This routine in this module is based on Bruno Sicardy’s program
   positionv.f. It gives the projected position (in the sky plane) of
   the instants for immersion and emersion for each observer. This
   routine depends on the output from previous modules and also has a
   strong dependence of a calculation for the sidereal time propagated
   to the occultation day for each site.

   It also depends on the geocentric position for each site, which in
   turn, depends on the geoid model adopted. Note that the model adopted
   is the WGS84, used by the GPS devices, which defines the position on
   the Earht surface on ITRS. In SORA we use astropy, which already
   takes into account the ITRS.

Ellipse Fit [SubSec:ellipse_fit]
--------------------------------

This module is an automated and updated module based on program
ellipse_fit.f, which is a program based on *Numerical Recipes* routines
to fit an ellipse with its 5 parameters on the data points (the
:math:`(\xi , \eta)_i` from the sky projection obtained from
positionv.f).

With the projections in the sky plane (:math:`\xi , \eta`) for each
immersion and emersion time for each station (observer and light curve)
we can find the best apparent ellipse. To describe an ellipse we need
five parameters: (i) and (ii) the ellipse centre (:math:`f_0,g_0`);
(iii) the apparent semi-major axis (:math:`a'`); (iv) the apparent
oblatness :math:`\left(\epsilon' = \dfrac{a' - b'}{a'}\right)`; and (v)
the position angle of the semi-minor axis (:math:`P`).

With less than 5 points (2.5 chords), it is not possible to fit all the
parameters to find the best apparent ellipse. However, some
considerations can help us. For instance, if some information about the
object is known (equivalent radius from thermal measures or previous
occultations, consider the body as spherical, etc), one or more
parameters can be constrained.

With one chord is only possible to determine the ellipse centre if the
radius is known. Since solutions are degenerate (usually called “North”
and “South” solutions), the result will be the two solutions plus an
average solution. If nothing is known the chord can indicate a limit to
the object’s size.

For version v0.1 the user have to provide the parameters to the fit
(initial guesses and ranges), since no automatic procedure was
implemented. After providing the points, SORA will fit an ellipse and
compare the radial distance of each point to the ellipse using a
:math:`\chi^2` method. This procedure is repeated for n-interactions
with the 5 parameters correlated and the minimum :math:`\chi^2` will
determine the best solution.

.. _SubSec:major_upgrades:

Major features to be implemented on future upgrades
---------------------------------------------------

Here we briefly present some of the intended modules to be implemented
in future versions of SORA: “body”, “rotation”, and “shape”.

The “body” will deal with the object characteristics. The “rotation”
will include features to work with rotational light curve and obtain
rotational periods. The “shape” is intended to work with how to obtain
the 3D shape (and its characteristics such as density) of an object
combining results from light curves, rotational period and occultation
data.

A module to deal with occultations by objects with atmosphere is also
planned, as well as the implementation of shortcuts for faster use and
the improvement of :math:`\chi^2` functions and parallelisation of the
procedures.

Note that SORA is an open source code so any new feature can be easily
implemented and tested and there is room for various implementations.

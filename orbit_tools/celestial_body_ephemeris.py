# Various Celestial Bodies in Earth Centered Inertial (ECI)
# Ryo Takatori
# 10/17/2020

# Modules
import numpy as np
import planet_data as pd

earth_radius = pd.Earth.radius
# Vallado Sun (Credits to David Vallado)
# %
# %------------------------------------------------------------------------------
# %
# %                           function sun
# %
# %  this function calculates the geocentric equatorial position vector
# %    the sun given the julian date.  this is the low precision formula and
# %    is valid for years from 1950 to 2050.  accuaracy of apparent coordinates
# %    is 0.01  degrees.  notice many of the calculations are performed in
# %    degrees, and are not changed until later.  this is due to the fact that
# %    the almanac uses degrees exclusively in their formulations.
# %
# %  author        : david vallado                  719-573-2600   27 may 2002
# %
# %  revisions
# %    vallado     - fix mean lon of sun                            7 mat 2004
# %
# %  inputs          description                    range / units
# %    jd          - julian date                    days from 4713 bc
# %
# %  outputs       :
# %    rsun        - ijk position vector of the sun au
# %    rtasc       - right ascension                rad
# %    decl        - declination                    rad
# %
# %  locals        :
# %    meanlong    - mean longitude
# %    meananomaly - mean anomaly
# %    eclplong    - ecliptic longitude
# %    obliquity   - mean obliquity of the ecliptic
# %    tut1        - julian centuries of ut1 from
# %                  jan 1, 2000 12h
# %    ttdb        - julian centuries of tdb from
# %                  jan 1, 2000 12h
# %    hr          - hours                          0 .. 24              10
# %    min         - minutes                        0 .. 59              15
# %    sec         - seconds                        0.0  .. 59.99          30.00
# %    temp        - temporary variable
# %    deg         - degrees
# %
# %  coupling      :
# %    none.
# %
# %  references    :
# %    vallado       2007, 281, alg 29, ex 5-1
# %
# % [rsun,rtasc,decl] = sun ( jd );
# % ------------------------------------------------------------------------------

def sun(jd):
    twopi = 2.0*np.pi

    # -------------------------  implementation   -----------------
    # -------------------  initialize values   --------------------
    tut1 = (jd-2451545.0)/ 36525.0

    meanlong = 280.460 + 36000.77*tut1
    meanlong = np.fmod(meanlong, 360.0) #deg

    ttdb = tut1
    meananomaly = 357.5277233 + 35999.05034 *ttdb
    meananomaly = np.fmod(np.deg2rad(meananomaly), twopi) #rad
    if meananomaly < 0.0:
        meananomaly = twopi + meananomaly

    eclplong = meanlong + 1.914666471*np.sin(meananomaly) \
                + 0.019994643 *np.sin(2.0*meananomaly) #deg
    eclplong = np.fmod(eclplong, 360.0) #deg

    obliquity = 23.439291 - 0.0130042 *ttdb #deg

    eclplong = np.deg2rad(eclplong)
    obliquity = np.deg2rad(obliquity)

    # --------- find magnitude of sun vector, )   components ------
    magr = 1.000140612  - 0.016708617*np.cos(meananomaly) \
                            - 0.000139589*np.cos(2.0*meananomaly)    # in au's

    rsun = np.empty(3)
    rsun[0] = magr*np.cos(eclplong)
    rsun[1] = magr*np.cos(obliquity)*np.sin(eclplong)
    rsun[2] = magr*np.sin(obliquity)*np.sin(eclplong)

    rtasc= np.arctan(np.cos(obliquity)*np.tan(eclplong))

    # --- check that rtasc is in the same quadrant as eclplong ----
    if eclplong < 0.0:
        eclplong = eclplong + twopi    # make sure it's in 0 to 2pi range
    if np.absolute(eclplong-rtasc) > np.pi*0.5:
        rtasc= rtasc + 0.5*np.pi*np.round((eclplong-rtasc)/(0.5*np.pi))
    decl = np.arcsin(np.sin(obliquity)*np.sin(eclplong))
    return rsun,rtasc,decl

# Curtis Moon (Credits to Curtis Howard)
def moon(jd):
    #
    #...Calculates the geocentric equatorial position vector of the moon
    # given the Julian day.
    #
    # User M-functions required: None
    # -------------------------------------------------------------------------

    # ------------------------- implementation -----------------
    #...Time in centuries since J2000:
    T = (jd-2451545)/36525
    #...Ecliptic longitude (deg):
    e_long = 218.32 + 481267.881*T \
    + 6.29*np.sin(np.deg2rad(135.0+477198.87*T)) - 1.27*np.sin(np.deg2rad(259.3-413335.36*T)) \
    + 0.66*np.sin(np.deg2rad(235.7+890534.22*T)) + 0.21*np.sin(np.deg2rad(269.9+954397.74*T)) \
    - 0.19*np.sin(np.deg2rad(357.5+35999.05*T)) - 0.11*np.sin(np.deg2rad(186.5+966404.03*T))
    e_long = np.mod(e_long, 360)
    e_long = np.deg2rad(e_long)
    #...Ecliptic latitude (deg):
    e_lat = 5.13*np.sin(np.deg2rad(93.3+483202.02*T)) + 0.28*np.sin(np.deg2rad(228.2+960400.89*T)) \
    - 0.28*np.sin(np.deg2rad(318.3+6003.15*T)) - 0.17*np.sin(np.deg2rad(217.6-407332.21*T))
    e_lat = np.mod(e_lat, 360)
    e_lat = np.deg2rad(e_lat)
    #...Horizontal parallax (deg):
    h_par = 0.9508 \
    + 0.0518*np.cos(np.deg2rad(135.0+477198.87*T)) + 0.0095*np.cos(np.deg2rad(259.3-413335.36*T)) \
    + 0.0078*np.cos(np.deg2rad(235.7+890534.22*T)) + 0.0028*np.cos(np.deg2rad(269.9+954397.74*T))
    h_par = np.mod(h_par,360)
    h_par = np.deg2rad(h_par)
    #...Angle between earth�s orbit and its equator (deg):
    obliquity = 23.439291 - 0.0130042*T
    obliquity = np.deg2rad(obliquity)
    #...Direction cosines of the moon�s geocentric equatorial position vector:
    l = np.cos(e_lat)*np.cos(e_long)
    m = np.cos(obliquity)*np.cos(e_lat)*np.sin(e_long) - np.sin(obliquity)*np.sin(e_lat)
    n = np.sin(obliquity)*np.cos(e_lat)*np.sin(e_long) + np.cos(obliquity)*np.sin(e_lat)
    #...Earth-moon distance (km):
    dist = earth_radius/np.sin(h_par)
    #...Moon�s geocentric equatorial position vector (km):
    r_moon = dist*np.array([l, m, n])
    return r_moon


# Planetary Ephemeris (Credits to Dr.A)



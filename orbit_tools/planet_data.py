# Planet Specifications
# Ryo Takatori
# 11/18/2020

# Modules
import numpy as np
from . import misc_data

au = misc_data.au
jc = misc_data.jc
j2000 = misc_data.j2000

# Data based from Curtis Textbook

class Sun:
    name = "Sun"
    radius = 696000 # radius [km]
    mass = 1.989e30 # mass [km]
    sidereal_rotation_period = 25.38*24 # sidereal rotation period [hour]
    inclination_to_equator = 7.25 # inclination of equator to orbit plane [deg]
    a = None # semi-major axis of orbit around Sun
    ecc = None # eccentricity of orbit around Sun
    inclination_to_ecliptic = None # inclination of orbit to the ecliptic plane [deg]
    orbit_sidereal_period = None # orbit sidereal period [day]
    mu = 132712000000 # gravitational parameter [km^3/s^2]
    soi_radius = None # sphere of influence radius
    mean_angular_rotation = 2*np.pi/sidereal_rotation_period/60/60 # mean angular rotation [rad/s]
    oblateness = None

    def eci_location(self, jd):
        # Vallado Sun (Credits to David Vallado)
        # #
        # #------------------------------------------------------------------------------
        # #
        # #                           function sun
        # #
        # #  this function calculates the geocentric equatorial position vector
        # #    the sun given the julian date.  this is the low precision formula and
        # #    is valid for years from 1950 to 2050.  accuaracy of apparent coordinates
        # #    is 0.01  degrees.  notice many of the calculations are performed in
        # #    degrees, and are not changed until later.  this is due to the fact that
        # #    the almanac uses degrees exclusively in their formulations.
        # #
        # #  author        : david vallado                  719-573-2600   27 may 2002
        # #
        # #  revisions
        # #    vallado     - fix mean lon of sun                            7 mat 2004
        # #
        # #  inputs          description                    range / units
        # #    jd          - julian date                    days from 4713 bc
        # #
        # #  outputs       :
        # #    rsun        - ijk position vector of the sun au
        # #    rtasc       - right ascension                rad
        # #    decl        - declination                    rad
        # #
        # #  locals        :
        # #    meanlong    - mean longitude
        # #    meananomaly - mean anomaly
        # #    eclplong    - ecliptic longitude
        # #    obliquity   - mean obliquity of the ecliptic
        # #    tut1        - julian centuries of ut1 from
        # #                  jan 1, 2000 12h
        # #    ttdb        - julian centuries of tdb from
        # #                  jan 1, 2000 12h
        # #    hr          - hours                          0 .. 24              10
        # #    min         - minutes                        0 .. 59              15
        # #    sec         - seconds                        0.0  .. 59.99          30.00
        # #    temp        - temporary variable
        # #    deg         - degrees
        # #
        # #  coupling      :
        # #    none.
        # #
        # #  references    :
        # #    vallado       2007, 281, alg 29, ex 5-1
        # #
        # # [rsun,rtasc,decl] = sun ( jd )
        # # ------------------------------------------------------------------------------
        twopi = 2.0*np.pi

        # -------------------------  implementation   -----------------
        # -------------------  initialize values   --------------------
        tut1 = (jd-2451545.0)/ 36525.0

        meanlong = 280.460 + 36000.77*tut1
        meanlong = np.fmod(meanlong, 360.0) #deg

        ttdb = tut1
        meananomaly = 357.5277233 + 35999.05034*ttdb
        meananomaly = np.fmod(np.deg2rad(meananomaly), twopi) #rad
        if meananomaly < 0.0:
            meananomaly = twopi + meananomaly

        eclplong = meanlong + 1.914666471*np.sin(meananomaly) \
                    + 0.019994643 *np.sin(2.0*meananomaly) #deg
        eclplong = np.fmod(eclplong, 360.0) #deg

        obliquity = 23.439291 - 0.0130042*ttdb #deg

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
            rtasc = rtasc + 0.5*np.pi*np.round((eclplong-rtasc)/(0.5*np.pi))
        decl = np.arcsin(np.sin(obliquity)*np.sin(eclplong))
        return rsun, rtasc, decl

class Mercury:
    name = "Mercury"
    radius = 2440 # radius [km]
    mass = 330.2e21 # mass [km]
    sidereal_rotation_period = 58.65*24 # sidereal rotation period [hour]
    inclination_to_equator = 0.01 # inclination of equator to orbit plane [deg]
    a = 57.91e6 # semi-major axis of orbit around Sun
    ecc = 0.2056 # eccentricity of orbit around Sun
    inclination_to_ecliptic = 7.00 # inclination of orbit to the ecliptic plane [deg]
    orbit_sidereal_period = 87.97 # orbit sidereal period [day]
    mu = 22030 # gravitational parameter [km^3/s^2]
    soi_radius = 112000 # sphere of influence radius
    mean_angular_rotation = 2*np.pi/sidereal_rotation_period/60/60 # mean angular rotation [rad/s]
    oblateness = 0.0000
    
    def ephemeris(self, t):
        t = (t-j2000)/jc
        a = 0.387098310 # AU but in km later
        ecc = 0.20563175 + 0.000020406*t - 0.0000000284*t**2 - 0.00000000017*t**3
        inc = 7.004986 - 0.0059516*t + 0.00000081*t**2 + 0.000000041*t**3 # degs
        raan = 48.330893 - 0.1254229*t - 0.00008833*t**2 - 0.000000196*t**3 # degs
        w_hat = 77.456119 + 0.1588643*t - 0.00001343*t**2 + 0.000000039*t**3 # degs
        L = 252.250906 + 149472.6746358*t - 0.00000535*t**2 + 0.000000002*t**3 # degs
        planet_coe = [a, ecc, inc, raan, w_hat, L]
        return planet_coe

class Venus:
    name = "Venus"
    radius = 6052 # radius [km]
    mass = 4.869e24 # mass [km]
    sidereal_rotation_period = 243*24 # sidereal rotation period [hour]
    inclination_to_equator = 177.4 # inclination of equator to orbit plane [deg]
    a = 108.2e6 # semi-major axis of orbit around Sun
    ecc = 0.0067 # eccentricity of orbit around Sun
    inclination_to_ecliptic = 3.39 # inclination of orbit to the ecliptic plane [deg]
    orbit_sidereal_period = 224.7 # orbit sidereal period [day]
    mu = 324900 # gravitational parameter [km^3/s^2]
    soi_radius = 616000 # sphere of influence radius
    mean_angular_rotation = 2*np.pi/sidereal_rotation_period/60/60 # mean angular rotation [rad/s]
    oblateness = 0.0000

    def ephemeris(self, t):
        t = (t-j2000)/jc
        a = 0.723329820 # AU
        ecc = 0.00677188 - 0.000047766*t + 0.000000097*t**2 + 0.00000000044*t**3
        inc = 3.394662 - 0.0008568*t - 0.00003244*t**2 + 0.000000010*t**3 # degs
        raan = 76.679920 - 0.2780080*t - 0.00014256*t**2 - 0.000000198*t**3 # degs
        w_hat = 131.563707 + 0.0048646*t - 0.00138232*t**2 - 0.000005332*t**3 # degs
        L = 181.979801 + 58517.8156760*t + 0.00000165*t**2 - 0.000000002*t**3 # degs
        planet_coe = [a, ecc, inc, raan, w_hat, L]
        return planet_coe

class Earth:
    name = "Earth"
    radius = 6378 # radius [km]
    mass = 5.974e24 # mass [km]
    sidereal_rotation_period = 23.9345 # sidereal rotation period [hour]
    inclination_to_equator = 23.45 # inclination of equator to orbit plane [deg]
    a = 149.6e6 # semi-major axis of orbit around Sun
    ecc = 0.0167 # eccentricity of orbit around Sun
    inclination_to_ecliptic = 0 # inclination of orbit to the ecliptic plane [deg]
    orbit_sidereal_period = 365.256 # orbit sidereal period [day]
    mu = 398600 # gravitational parameter [km^3/s^2]
    soi_radius = 925000 # sphere of influence radius
    mean_angular_rotation = 2*np.pi/sidereal_rotation_period/60/60 # mean angular rotation [rad/s]
    oblateness = 0.003353

    def ephemeris(self, t):
        t = (t-j2000)/jc
        a = 1.000001018 # AU
        ecc = 0.01670862 - 0.000042037*t - 0.0000001236*t**2 + 0.00000000004*t**3
        inc = 0.0000000 + 0.0130546*t - 0.00000931*t**2 - 0.000000034*t**3 # degs
        raan = 0.0 # degs
        w_hat = 102.937348 + 0.3225557*t + 0.00015026*t**2 + 0.000000478*t**3 # degs
        L = 100.466449 + 35999.372851*t - 0.00000568*t**2 + 0.000000000*t**3 # degs
        planet_coe = [a, ecc, inc, raan, w_hat, L]
        return planet_coe

class Moon:
    name = "Moon"
    radius = 1737 # radius [km]
    mass = 73.48e21 # mass [km]
    sidereal_rotation_period = 27.32*24 # sidereal rotation period [hour]
    inclination_to_equator = 6.68 # inclination of equator to orbit plane [deg]
    a = 384.4e3 # semi-major axis of orbit around Sun
    ecc = 0.0549 # eccentricity of orbit around Sun
    inclination_to_ecliptic = 5.145 # inclination of orbit to the ecliptic plane [deg]
    orbit_sidereal_period = 27.322 # orbit sidereal period [day]
    mu = 4903 # gravitational parameter [km^3/s^2]
    soi_radius = 66100 # sphere of influence radius
    mean_angular_rotation = 2*np.pi/sidereal_rotation_period/60/60 # mean angular rotation [rad/s]
    oblateness = 0.0012

    # Curtis Moon (Credits to Curtis Howard)
    def eci_location(self, jd):
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
        dist = Earth().radius/np.sin(h_par)
        #...Moon's geocentric equatorial position vector (km):
        r_moon = dist*np.array([l, m, n])
        return r_moon

class Mars:
    name = "Mars"
    radius = 3396 # radius [km]
    mass = 641.9e21 # mass [km]
    sidereal_rotation_period = 24.62 # sidereal rotation period [hour]
    inclination_to_equator = 25.19 # inclination of equator to orbit plane [deg]
    a = 227.9e6 # semi-major axis of orbit around Sun
    ecc = 0.0935 # eccentricity of orbit around Sun
    inclination_to_ecliptic = 1.850 # inclination of orbit to the ecliptic plane [deg]
    orbit_sidereal_period = 1.881*365 # orbit sidereal period [day]
    mu = 42828 # gravitational parameter [km^3/s^2]
    soi_radius = 577000 # sphere of influence radius
    mean_angular_rotation = 2*np.pi/sidereal_rotation_period/60/60 # mean angular rotation [rad/s]
    oblateness = 0.00648

    def ephemeris(self, t):
        t = (t-j2000)/jc
        a = 1.523679342 # AU
        ecc = 0.09340062 + 0.000090483*t - 0.00000000806*t**2 - 0.00000000035*t**3
        inc = 1.849726 - 0.0081479*t - 0.00002255*t**2 - 0.000000027*t**3 # degs
        raan = 49.558093 - 0.2949846*t - 0.00063993*t**2 - 0.000002143*t**3 # degs
        w_hat = 336.060234 + 0.4438898*t - 0.00017321*t**2 + 0.000000300*t**3 # degs
        L = 355.433275 + 19140.2993313*t + 0.00000261*t**2 - 0.000000003*t**3 # degs
        planet_coe = [a, ecc, inc, raan, w_hat, L]
        return planet_coe

class Jupiter:
    name = "Jupiter"
    radius = 71490 # radius [km]
    mass = 1.899e27 # mass [km]
    sidereal_rotation_period = 9.925 # sidereal rotation period [hour]
    inclination_to_equator = 3.13 # inclination of equator to orbit plane [deg]
    a = 778.6e6 # semi-major axis of orbit around Sun
    ecc = 0.0489 # eccentricity of orbit around Sun
    inclination_to_ecliptic = 1.304 # inclination of orbit to the ecliptic plane [deg]
    orbit_sidereal_period = 11.86*365 # orbit sidereal period [day]
    mu = 126686000 # gravitational parameter [km^3/s^2]
    soi_radius = 48200000 # sphere of influence radius
    mean_angular_rotation = 2*np.pi/sidereal_rotation_period/60/60 # mean angular rotation [rad/s]
    oblateness = 0.06487

    def ephemeris(self, t):
        t = (t-j2000)/jc
        a = 5.202603191 + 0.0000001913*t # AU
        ecc = 0.04849485 + 0.000163244*t - 0.0000004719*t**2 + 0.00000000197*t**3
        inc = 1.303270 - 0.0019872*t + 0.00003318*t**2 + 0.000000092*t**3 # degs
        raan = 100.464441 + 0.1766828*t + 0.00090387*t**2 - 0.000007032*t**3 # degs
        w_hat = 14.331309 + 0.2155525*t + 0.00072252*t**2 - 0.000004590*t**3 # degs
        L = 34.351484 + 3034.9056746*t - 0.00008501*t**2 + 0.000000004*t**3 # degs
        planet_coe = [a, ecc, inc, raan, w_hat, L]
        return planet_coe

class Saturn:
    name = "Saturn"
    radius = 60270 # radius [km]
    mass = 568.5e4 # mass [km]
    sidereal_rotation_period = 10.66 # sidereal rotation period [hour]
    inclination_to_equator = 26.73 # inclination of equator to orbit plane [deg]
    a = 1.433e9 # semi-major axis of orbit around Sun
    ecc = 0.0565 # eccentricity of orbit around Sun
    inclination_to_ecliptic = 2.485 # inclination of orbit to the ecliptic plane [deg]
    orbit_sidereal_period = 29.46*365 # orbit sidereal period [day]
    mu = 37931000 # gravitational parameter [km^3/s^2]
    soi_radius = 54800000 # sphere of influence radius
    mean_angular_rotation = 2*np.pi/sidereal_rotation_period/60/60 # mean angular rotation [rad/s]
    oblateness = 0.09796

    def ephemeris(self, t):
        t = (t-j2000)/jc
        a = 9.5549009596 - 0.0000021389*t # AU
        ecc = 0.05550862 - 0.000346818*t - 0.0000006456*t**2 + 0.00000000338*t**3
        inc = 2.488878 + 0.0025515*t - 0.00004903*t**2 + 0.000000018*t**3 # degs
        raan = 113.665524 - 0.2566649*t - 0.00018345*t**2 + 0.000000357*t**3 # degs
        w_hat = 93.056787 + 0.5665496*t + 0.00052809*t**2 - 0.000004882*t**3 # degs
        L = 50.077471 + 1222.1137943*t + 0.00021004*t**2 - 0.000000019*t**3 # degs
        planet_coe = [a, ecc, inc, raan, w_hat, L]
        return planet_coe

class Uranus:
    name = "Uranus"
    radius = 25560 # radius [km]
    mass = 86.83e24 # mass [km]
    sidereal_rotation_period = 17.24 # sidereal rotation period [hour]
    inclination_to_equator = 97.77 # inclination of equator to orbit plane [deg]
    a = 2.872e9 # semi-major axis of orbit around Sun
    ecc = 0.0457 # eccentricity of orbit around Sun
    inclination_to_ecliptic = 0.772 # inclination of orbit to the ecliptic plane [deg]
    orbit_sidereal_period = 84.01*365 # orbit sidereal period [day]
    mu = 5794000 # gravitational parameter [km^3/s^2]
    soi_radius = 51800000 # sphere of influence radius
    mean_angular_rotation = 2*np.pi/sidereal_rotation_period/60/60 # mean angular rotation [rad/s]
    oblateness = 0.02293

    def ephemeris(self, t):
        t = (t-j2000)/jc
        a = 19.218446062 - 0.0000000372*t + 0.00000000098*t**2 # AU
        ecc = 0.04629590 - 0.000027337*t + 0.0000000790*t**2 + 0.00000000025*t**3
        inc = 0.773196 - 0.0016869*t + 0.00000349*t**2 + 0.00000000016*t**3 # degs
        raan = 74.005947 + 0.0741461*t + 0.00040540*t**2 + 0.000000104*t**3 # degs
        w_hat = 173.005159 + 0.0893206*t - 0.00009470*t**2 + 0.000000413*t**3 # degs
        L = 314.055005 + 428.4669983*t - 0.00000486*t**2 - 0.000000006*t**3 # degs
        planet_coe = [a, ecc, inc, raan, w_hat, L]
        return planet_coe

class Neptune:
    name = "Neptune"
    radius = 4760 # radius [km]
    mass = 102.4e24 # mass [km]
    sidereal_rotation_period = 16.11 # sidereal rotation period [hour]
    inclination_to_equator = 28.32 # inclination of equator to orbit plane [deg]
    a = 4.495e9 # semi-major axis of orbit around Sun
    ecc = 0.0113 # eccentricity of orbit around Sun
    inclination_to_ecliptic = 1.769 # inclination of orbit to the ecliptic plane [deg]
    orbit_sidereal_period = 164.8*365 # orbit sidereal period [day]
    mu = 6835100 # gravitational parameter [km^3/s^2]
    soi_radius = 86600000 # sphere of influence radius
    mean_angular_rotation = 2*np.pi/sidereal_rotation_period/60/60 # mean angular rotation [rad/s]
    oblateness = 0.01708

    def ephemeris(self, t):
        t = (t-j2000)/jc
        a = 30.110386869 - 0.0000001663*t + 0.00000000069*t**2 # AU
        ecc = 0.00898809 + 0.000006408*t - 0.0000000008*t**2
        inc = 1.769952 + 0.0002557*t + 0.00000023*t**2 - 0.0000000000*t**3 # degs
        raan = 131.784057 - 0.0061651*t - 0.00000219*t**2 - 0.000000078*t**3 # degs
        w_hat = 48.123691 + 0.0291587*t + 0.00007051*t**2 - 0.000000000*t**3 # degs
        L = 304.348665 + 218.4862002*t + 0.00000059*t**2 - 0.000000002*t**3 # degs
        planet_coe = [a, ecc, inc, raan, w_hat, L]
        return planet_coe

class Pluto:
    name = "Pluto"
    radius = 1195 # radius [km]
    mass = 12.5e21 # mass [km]
    sidereal_rotation_period = 6.387*24 # sidereal rotation period [hour]
    inclination_to_equator = 122.5 # inclination of equator to orbit plane [deg]
    a = 5.870e9 # semi-major axis of orbit around Sun
    ecc = 0.2444 # eccentricity of orbit around Sun
    inclination_to_ecliptic = 17.16 # inclination of orbit to the ecliptic plane [deg]
    orbit_sidereal_period = 247.7*365 # orbit sidereal period [day]
    mu = 830 # gravitational parameter [km^3/s^2]
    soi_radius = 3080000 # sphere of influence radius
    mean_angular_rotation = 2*np.pi/sidereal_rotation_period/60/60 # mean angular rotation [rad/s]
    oblateness = None
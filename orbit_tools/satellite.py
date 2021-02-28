# Class for Satellite
# Ryo Takatori
# 11/17/2020

# Modules
import numpy as np
import datetime
import tle_reader as tr
import coe
import planet_data as pd

earth_mu = pd.Earth.mu # Earth gravitational constant [km^3/s^2]
earth_radius = pd.Earth.radius # Earth radius [km]

class Satellite:
    def __init__(self, filename):
        # Two-Line Element Set
        self.tle = tr.tle_reader(filename)
        self.name = self.tle[0][0] # name
        # Line 1
        self.line_number_1 = self.tle[1][0] # line number
        self.catalog_number_1 = self.tle[1][1] # satellite catalog number
        self.classification = self.tle[1][2] # classification
        # International designator
        self.id_launch_year = self.tle[1][3] # last two digits of launch year
        self.id_launch_number = self.tle[1][4] # launch number of the year
        self.id_launch_piece = self.tle[1][5] # piece of launch
        #
        self.epoch_year = self.tle[1][6] # epoch year (last two digits of launch year)
        self.epoch_day_of_year = self.tle[1][7] # epoch (day of the year and fractional portion of the day)
        self.ballistic_coefficient = self.tle[1][8] # first derivative of mean motion (ballistic coefficient)
        self.d2_mean_motion = self.tle[1][9] # second derivative of mean motion
        self.bstar = self.tle[1][10] # BSTAR
        self.ephemeris_type = self.tle[1][11] # ephemeris type
        self.element_set_number = self.tle[1][12] # element set number
        self.checksum_1 = self.tle[1][13] # check sum
        # Line 2
        self.line_number_2 = self.tle[2][0] # line number
        self.catalog_number_2 = self.tle[2][1] # satellite catalog number
        self.inclination = self.tle[2][2] # inclination [deg]
        self.right_ascension_ascending_node = self.tle[2][3] # right ascension of ascending node [deg]
        self.eccentricity = self.tle[2][4] # eccentricity
        self.argument_of_perigee = self.tle[2][5] # argument of perigee [deg]
        self.mean_anomaly = self.tle[2][6] # mean anomaly [deg]
        self.mean_motion = self.tle[2][7] # mean motion [rev/day]
        self.revolution = self.tle[2][8] # revolution number at epoch [rev]
        self.checksum_2 = self.tle[2][9] # check sum

        def ecc_anomaly_newtons_method(ecc, me):
            if me < np.pi:
                E = me + (ecc/2)
            elif me > np.pi:
                E = me - (ecc/2)
            else:
                pass
            x = 1
            tol = 1
            while tol > 10e-8:
                if x < 1000:
                    E = E + ((me-E+ecc*np.sin(E))/(1-ecc*np.cos(E)))
                tol = ((me-E+ecc*np.sin(E))/(1-ecc*np.cos(E)))
                x += 1
            return E
        # Classical Orbital Elements and other related elements and parameters
        self.n = (self.mean_motion*2*np.pi)/(24*60*60) # mean motion
        self.me = self.mean_anomaly # mean anomaly [deg]
        self.a = (earth_mu/self.n**2)**(1/3) # semi-major axis [km]
        self.T = 2*np.pi*np.sqrt((self.a**3)/earth_mu) # period [s]
        self.ecc = self.eccentricity # eccentricity
        self.rp = self.a*(1-self.ecc) # perigee [km]
        self.ra = (2*self.a) - self.rp # apogee [km]
        self.zp = self.rp - earth_radius # perigee altitude [km]
        self.za = self.ra - earth_radius # apogee altitude [km]
        self.E = ecc_anomaly_newtons_method(self.ecc, self.me) # eccentric anomaly
        self.h = np.sqrt(self.a*earth_mu*(1-self.ecc**2)) # specifc engular momentum [km^2/s]
        self.inc = self.inclination # inclination [deg]
        self.raan = self.right_ascension_ascending_node # right ascension of ascending node [deg]
        self.w = self.argument_of_perigee # argument of perigee [deg]
        self.theta = np.rad2deg(2*np.arctan(np.sqrt(((1+self.ecc)/(1-self.ecc)))*np.tan(self.E/2))) # true anomaly [deg]
        self.coe = [self.a, self.ecc, self.h, self.inc, self.raan, self.w, self.theta]
        # Vector
        r, v = coe.coe2rv(earth_mu, self.coe, False) # coe to vectors
        self.r = r # position vector [km]
        self.v = v # velocity vector [km/s]
        # Check for past 2000 or before
        if self.epoch_year > 56:
            epoch_full_year = 1900 + self.epoch_year
        else:
            epoch_full_year = 2000 + self.epoch_year
        self.epoch = datetime.datetime(epoch_full_year, 1, 1) + datetime.timedelta(self.epoch_day_of_year - 1) # epoch in standard datetime
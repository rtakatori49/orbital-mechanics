import datetime

class MiscellaneousData:
    def __init__(self):
        self.c = 2.998e8 # speed of light [m/s]
        self.au = 149597870 # 1 astronomical unit to kilometer
        self.j4_doy = datetime.datetime(2020, 7, 4, 0, 0, 0).timetuple().tm_yday # July 4th in day of year
        self.j2000 = 2451545 # J2000 in julian date
        self.jc = 36525 # Julian century
        self.space_alt = 100 # apace altitude [km]
        self.c_d = 2.2 # commmon coefficient of drag for spacecraft
        self.c_r = 1.2 # common coefficient of reflection for spacecraft
misc_data = MiscellaneousData()
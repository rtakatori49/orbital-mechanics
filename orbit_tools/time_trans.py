# Time Transformations
# Ryo Takatori
# 10/16/2020

# Modules
import datetime
import numpy as np
from . import planet_data as p_d

# Julian Date
def jd(input_datetime):
    # Extract date information
    y = input_datetime.year # year
    m = input_datetime.month # month
    d = input_datetime.day # day
    hr = input_datetime.hour # hour
    mn = input_datetime.minute # minute
    s = input_datetime.second # second

    # Julian Date
    j0 = 367*y \
        - np.floor((7*(y+np.floor((m + 9)/12)))/4) \
        + np.floor((275*m)/9) + d + 1721013.5
    jd = j0 + (hr/24) + (mn/(24*60)) + (s/(24*60*60)) # add time
    return jd

# Local Sidereal Time
def lst(input_datetime, lon):
    # Make angle between 0 and 360
    def deg_mod(angle):
        if angle < 0:
            angle %= -360
        else:
            angle %= 360
        return angle
    # Extract date information
    y = input_datetime.year # year
    m = input_datetime.month # month
    d = input_datetime.day # day
    hr = input_datetime.hour # hour
    mn = input_datetime.minute # minute
    s = input_datetime.second # second
    input_date = datetime.datetime(y, m, d, 0, 0, 0)
    ut1 = [hr, mn, s]
    w_earth = np.rad2deg(p_d.Earth().mean_angular_rotation)
    # Date calculation
    jd0 = jd(input_date)
    # Time calculation
    t_ut1 = (jd0-2451545.0) / 36525
    # Greenwich sidereal time at 0hr
    theta_gst0 = 100.4606184 + 36000.77004*t_ut1 + 0.000387933*(t_ut1)**2 \
        - 2.583e-8*(t_ut1)**3
    theta_gst0 = deg_mod(theta_gst0) # make GST0 between 0 and 360
    theta_gst = theta_gst0 + w_earth*(ut1[0]*60*60+ut1[1]*60+ut1[2])
    theta_lst = theta_gst + lon # local sidereal time
    theta_lst = deg_mod(theta_lst) # make LST between 0 and 360
    return theta_gst, theta_lst
# Time Transformations
# Ryo Takatori
# 10/16/2020

# Modules
import numpy as np
import datetime

# Julian Date
def jd(input_date):
    # Extract date information
    y = input_date.year # year
    m = input_date.month # month
    d = input_date.day # day
    hr = input_date.hour # hour
    mn = input_date.minute # minute
    s = input_date.second # second

    # Julian Date
    j_0 = 367*y - np.floor((7*(y+np.floor((m + 9)/12)))/4) + np.floor((275*m)/9) + d + 1721013.5
    jd = j_0 + (hr/24) + (mn/(24*60)) + (s/(24*60*60)) # add time
    return jd

# Local Sidereal Time
def lst():
    pass
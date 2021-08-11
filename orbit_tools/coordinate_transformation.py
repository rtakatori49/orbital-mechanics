# Coordinate Tranformation
# Ryo Takatori
# 10/16/2020

# Modules
import numpy as np
from . import direction_cosine_matrix as dcm
from . import planet_data as p_d
from . import time_trans as t_t

# Perifocal to Earth Centered Inertial (ECI)
def perifocal_to_eci(r_peri, v_peri, w, raan, inc):
    # Rotation 3 through argument of perigee
    cz_w = dcm.rot_z(w)
    cz_w[1, 0] = -cz_w[1, 0]
    cz_w[0, 1] = -cz_w[0, 1]

    # Rotation 1 through inclination
    cx_inc = dcm.rot_x(inc)
    cx_inc[1, 2] = -cx_inc[1, 2]
    cx_inc[2, 1] = -cx_inc[2, 1]

    # Rotation 3 through RAAN
    cz_raan = dcm.rot_z(raan)
    cz_raan[1, 0] = -cz_raan[1, 0]
    cz_raan[0, 1] = -cz_raan[0, 1]

    # DCM
    q = np.dot(cz_w, np.dot(cx_inc, cz_raan))
    q = q.T

    # States
    r = np.dot(q, r_peri.T)
    v = np.dot(q, v_peri.T)
    return r, v

# Earth Centered Inertial to Local Vertical Local Horizontal
# (ECI to LVLH)
def eci_to_lvlh(r, v):
    r = np.array(r)
    v = np.array(v)
    norm_r = np.linalg.norm(r) # Position magnitude [km]
    h = np.cross(r, v)
    norm_h = np.linalg.norm(h)
    i_hat = r/norm_r # i hat
    k_hat = h/norm_h # k hat
    j_hat = np.cross(k_hat, i_hat) # j hat
    q = np.array([i_hat, j_hat, k_hat]) # ECI to LVLH DCM
    r_lvlh = np.dot(q, r) # Relative position in LVLH [km]
    return r_lvlh, q

# Latitude Longitude Altitude to Earth Centered Inertial (LLA to ECI)
def lla_to_eci(lat, long, alt, date_time, geodetic=True):
    earth = p_d.Earth()
    f = earth.oblateness
    re = earth.radius
    lat = np.deg2rad(lat)
    _, theta_lst = t_t.lst(date_time, long)
    theta_lst = np.deg2rad(theta_lst)
    if geodetic:
        num = (1-f)**2
        den = np.sqrt(1-(2*f-(f**2)*np.sin(lat)**2))
        r = (re/den) + alt
        r_z = ((re*num)/den) + alt
        z = r_z*np.sin(lat)
    else:
        r = re + alt
        z = r*np.sin(lat)
    x = r*np.cos(lat)*np.cos(theta_lst)
    y = r*np.cos(lat)*np.sin(theta_lst)
    return np.array([x, y, z])

# Make angle between 0 and 360
def deg_mod(angle):
    if angle < 0:
        angle %= -360
    else:
        angle %= 360
    return angle

# Right Asecention and Declination to Azimuth and Elevation
def radec_to_azel(ra, dec, lat, long, date_time):
    dec = np.deg2rad(dec)
    lat = np.deg2rad(lat)
    long = np.deg2rad(long)
    _, theta_lst = t_t.lst(date_time, long)
    lha = np.deg2rad(deg_mod(theta_lst-ra))
    elevation = np.arcsin(np.sin(lat)*np.sin(dec)
        +np.cos(lat)*np.cos(dec)*np.cos(lha))
    azimuth = np.arctan2(-np.sin(lha)*np.cos(dec)/np.cos(elevation),
        (np.sin(dec)-np.sin(elevation)*np.sin(lat))/(np.cos(elevation)
        *np.cos(lat)))
    return np.rad2deg(azimuth), np.rad2deg(elevation)

def azel_to_radec(az, el, lat, long, date_time):
    az = np.deg2rad(az)
    el = np.deg2rad(el)
    lat = np.deg2rad(lat)
    _, theta_lst = t_t.lst(date_time, long)
    dec = np.arcsin(np.sin(el)*np.sin(lat)+np.cos(el)*np.cos(lat)*np.cos(az))
    lha = np.rad2deg(np.arctan2(-np.sin(az)*np.cos(el)/np.cos(dec),
    (np.sin(el)-np.sin(dec)*np.sin(lat))/(np.cos(dec)*np.cos(lat))))
    ra = deg_mod(theta_lst-lha)
    return ra, np.rad2deg(dec)

def state_to_razel(r, date_time, lat, long, alt):
    r_site = lla_to_eci(lat, long, alt, date_time)
    rho = r - r_site
    rho_norm = np.linalg.norm(rho)
    rho_hat = rho/rho_norm
    dec = np.rad2deg(np.arcsin(rho_hat[2]))
    ra = np.rad2deg(np.arccos(rho_hat[0]/np.cos(np.deg2rad(dec))))
    if rho[1] <= 0:
        ra = 360 - ra
    az, el = radec_to_azel(ra, dec, lat, long, date_time)
    return rho_norm, az, el
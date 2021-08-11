# Perturbations
# Ryo Takatori
# 10/16/2020

# Modules
import numpy as np
import datetime
from . import time_trans as t_t
from . import planet_data as p_d
from . import misc_data

# Constants
earth = p_d.Earth()
mue = earth.mu # Earth gravitational constant [km^3/s^2]
re = earth.radius # Earth radius [km]
earth_ang_vel = np.array([0, 0, earth.mean_angular_rotation]) # Earth Angular velocity [rad/s]
au = misc_data.au

# Class for perturbation definition before running equation of motion
class Perturbation:
    def __init__(self):
        self.drag = False
        self.srp = False
        self.n_body = False
        self.j2_6 = False
        self.area = 0
        self.mass = 0
        self.c_d = misc_data.c_d
        self.c_r = misc_data.c_r
        self.epoch = 0
        self.planets = []
        self.j = 6
        self.radius = re
        self.alt = misc_data.space_alt
        
    def set_drag(self, area, mass, c_d):
        self.drag = True
        self.area = area
        self.mass = mass
        self.c_d = c_d

    def set_srp(self, area, mass, c_r, epoch):
        self.srp = True
        self.area = area
        self.mass = mass
        self.c_d = c_r
        self.epoch = epoch

    def set_n_body(self, epoch, planets):
        self.n_body = True
        self.epoch = epoch
        self.planets = planets

    def set_j2_6(self, j):
        self.j2_6 = True
        self.j = j

    def set_all(self, area, mass, c_d, c_r, epoch, planets, j):
        self.set_drag(area, mass, c_d)
        self.set_srp(area, mass, c_r, epoch)
        self.set_n_body(epoch, planets)
        self.set_j2_6(j)

# Exponential Drag
def drag_exp(r, v, perts):
    A = perts.area
    m = perts.mass
    c_d = perts.c_d
    r = np.array(r)
    v = np.array(v)
    # Relative velocity [km/s]
    v_rel = v - np.cross(earth_ang_vel, r)
    norm_v_rel = np.linalg.norm(v_rel)
    # Altitude [km]
    norm_r = np.linalg.norm(r)
    alt = norm_r - re
    # Density [kg/m^3]
    rho = atmosphere(alt)
    a = -(0.5*c_d*A*rho*norm_v_rel*v_rel*1000)/m
    return a

# 1976 US Standard Atmosphere by Curtis
def atmosphere(z):
    # Geometric altitudes [km]
    h = [0, 25, 30, 40, 50, 60, 70, \
        80, 90, 100, 110, 120, 130, 140, \
        150, 180, 200, 250, 300, 350, 400, \
        450, 500, 600, 700, 800, 900, 1000]

    # Corresponding densities [kg/m^3] from USSA76
    r = [1.225, 4.008e-2, 1.841e-2, 3.996e-3, 1.027e-3, 3.097e-4, 8.283e-5, \
        1.846e-5, 3.416e-6, 5.606e-7, 9.708e-8, 2.222e-8, 8.152e-9, 3.831e-9, \
        2.076e-9, 5.194e-10, 2.541e-10, 6.073e-11, 1.916e-11, 7.014e-12, 2.803e-12, \
        1.184e-12, 5.215e-13, 1.137e-13, 3.070e-14, 1.136e-14, 5.759e-15, 3.561e-15]

    # Scale heights [km]
    H = [7.310, 6.427, 6.546, 7.360, 8.342, 7.583, 6.661, \
        5.927, 5.533, 5.703, 6.782, 9.973, 13.243, 16.322, \
        21.652, 27.974, 34.934, 43.342, 49.755, 54.513, 58.019, \
        60.980, 65.654, 76.377, 100.587, 147.203, 208.020]

    # Handle altitudes outside of the range
    if z > 1000:
        z = 1000
    elif z < 0:
        z = 0

    # Determine the interpolation interval:
    for j in range(27):
        if z >= h[j] and z < h[j+1]:
            i = j
    if z == 1000:
        i = 26
    # Exponential interpolation:
    density = r[i]*np.exp(-(z-h[i])/H[i])
    return density

# Solar Radiation Pressure
def srp(t, r, perts):
    A = perts.area
    m = perts.mass
    c_r = perts.c_r
    epoch = perts.epoch
    r = np.array(r)
    norm_r = np.linalg.norm(r)
    # Time adjustment
    current_date = epoch + datetime.timedelta(seconds=t) # current datetime
    doy = current_date.timetuple().tm_yday # current time in date of year
    dsa = doy - misc_data.j4_doy # days since aphelion (July 4)
    if dsa < 0:
        dsa = dsa + 365
    d_a = 2*np.pi*dsa
    s = 1358/(1.004 + 0.0534*np.cos(d_a)) # solar flux calculation [W/m^2]
    p_sr = s/misc_data.c
    jd = t_t.jd(epoch) # Julian date
    nu = shadow(jd, t, r) # check shadow
    a = -((p_sr*c_r*A)/m)*(r/norm_r)*nu
    return a

# Shadow function
def shadow(jd, t, r):
    sun = p_d.Sun()
    [r_sun, _, _] = sun.eci_location(jd+t/(24*60*60)) # Vallado sun calculation
    # [(sun vector) [km],(RAAN) [rad],(decl) [rad]]
    norm_r = np.linalg.norm(r) # Spacecraft position magnitude [km]
    norm_r_sun = np.linalg.norm(r_sun)*au # Sun position magnitude [km]
    # Angles from center of Earth required to calculate shadow
    theta = np.arccos(np.dot(r_sun*au,r)/(norm_r_sun*norm_r))
    theta_1 = np.arccos(re/norm_r)
    theta_2 = np.arccos(re/norm_r_sun)
    # Checking shadow
    if theta_1 + theta_2 <= theta:
        nu = 0
    else:
        nu = 1
    return nu

# N-body
def n_body(r, r_body, mu_body):
    r = np.array(r)
    norm_r_body = np.linalg.norm(r_body)
    # Satellite to body
    r_s_b = r_body - r
    norm_r_s_b = np.linalg.norm(r_s_b)

    # Binomial expansion
    q = np.dot(r, ((2*r_body) - r))/norm_r_body**2
    F = (((q**2)-(3*q)+3)/(1+((1-q)**(3/2))))*q

    # Perturbation
    a = (mu_body/(norm_r_s_b**3))*((F*r_body)-r)
    return a

# J2-6
def j2_6(r, perts):
    r = np.array(r)
    # Zonal harmonics
    j2 = 0.00108263
    j3 = -2.33936e-3*j2
    j4 = -1.49601e-3*j2
    j5 = -0.20995e-3*j2
    j6 = 0.49941e-3*j2

    norm_r = np.linalg.norm(r) # Position magnitude [km]

    # Perturbation calculation
    # J2
    a_J2 = np.empty(3)
    factor_J2 = -(3*j2*mue*re**2)/(2*norm_r**5)
    a_J2[0] = factor_J2*r[0]*(1 - ((5*(r[2]**2))/norm_r**2))
    a_J2[1] = factor_J2*r[1]*(1 - ((5*(r[2]**2))/norm_r**2))
    a_J2[2] = factor_J2*r[2]*(3 - ((5*(r[2]**2))/norm_r**2))

    # J3
    a_J3 = np.empty(3)
    factor_J3 = -(5*j3*mue*re**3)/(2*norm_r**7)
    a_J3[0] = factor_J3*r[0]*((3*r[2]) - ((7*(r[2]**3))/norm_r**2))
    a_J3[1] = factor_J3*r[1]*((3*r[2]) - ((7*(r[2]**3))/norm_r**2))
    a_J3[2] = factor_J3*(6*r[2]**2 - ((7*(r[2]**4))/norm_r**2) - ((3/5)*norm_r**2))

    # J4
    a_J4 = np.empty(3)
    factor_J4 = (15*j4*mue*re**4)/(8*norm_r**7)
    a_J4[0] = factor_J4*r[0]*(1 - ((14*(r[2]**2))/norm_r**2) + ((21*(r[2]**4))/norm_r**4))
    a_J4[1] = factor_J4*r[1]*(1 - ((14*(r[2]**2))/norm_r**2) + ((21*(r[2]**4))/norm_r**4))
    a_J4[2] = factor_J4*r[2]*(5 - ((70*(r[2]**2))/(3*norm_r**2)) + ((21*(r[2]**4))/norm_r**4))

    # J5
    a_J5 = np.empty(3)
    factor_J5 = (3*j5*mue*re**5*r[2])/(8*norm_r**9)
    a_J5[0] = factor_J5*r[0]*(35 - ((210*(r[2]**2))/norm_r**2) + ((231*(r[2]**4))/norm_r**4))
    a_J5[1] = factor_J5*r[1]*(35 - ((210*(r[2]**2))/norm_r**2) + ((231*(r[2]**4))/norm_r**4))
    a_J5[2] = factor_J5*r[2]*(105 - ((315*(r[2]**2))/norm_r**2) + ((231*(r[2]**4))/norm_r**4)) \
        - (15*j5*mue*re**5)/(8*norm_r**7)

    # J6
    a_J6 = np.empty(3)
    factor_J6 = -(j6*mue*re**6)/(16*norm_r**9)
    a_J6[0] = factor_J6*r[0]*(35 - ((945*(r[2]**2))/norm_r**2) + ((3465*(r[2]**4))/norm_r**4) \
        - ((3003*(r[2]**6))/norm_r**6))
    a_J6[1] = factor_J6*r[1]*(35 - ((945*(r[2]**2))/norm_r**2) + ((3465*(r[2]**4))/norm_r**4) \
        - ((3003*(r[2]**6))/norm_r**6))
    a_J6[2] = factor_J6*r[2]*(245 - ((2205*(r[2]**2))/(3*norm_r**2)) \
        + ((4851*(r[2]**4))/norm_r**4) - ((3003*(r[2]**6))/norm_r**6))

    # Perturbation
    a_J = [a_J2, a_J3, a_J4, a_J5, a_J6]
    a = sum(a_J[0:perts.j-1])
    #a = a_J2 + a_J3 + a_J4 + a_J5 + a_J6
    return a

def acceleration(t, r, v, perts):
    a_zero = np.array([0,0,0])
    # Drag
    if perts.drag:
        a_drag = drag_exp(r, v, perts)
    else:
        a_drag = a_zero
    # Solar radiation pressure
    if perts.srp:
        a_srp = srp(t, r, perts)
    else:
        a_srp = a_zero
    # N-body effects
    if perts.n_body:
        a_n_body_list = []
        # Account for all planets in list
        for planet in perts.planets:
            jd = t_t.jd(perts.epoch) # julian date
            # Moon
            if planet == "Moon":
                moon = p_d.Moon()
                r_body = moon.eci_location(jd+(t/(24*60*60)))
                mu_body = moon.mu
            # Sun
            if planet == "Sun":
                sun = p_d.Sun()
                [r_sun, _, _] = sun.eci_location(jd+(t/(24*60*60)))
                r_body = r_sun*au
                mu_body = sun.mu
            a_n_body_temp = n_body(r, r_body, mu_body)
            a_n_body_list.append(a_n_body_temp)
        a_n_body = sum(a_n_body_list)
    else:
        a_n_body = a_zero
    # J2~6
    if perts.j2_6:
        a_j2_6 = j2_6(r, perts)
    else:
        a_j2_6 = a_zero
    ax, ay, az = a_drag + a_srp + a_n_body + a_j2_6 # acceleration
    return [ax, ay, az]
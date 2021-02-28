# Classical Orbital Elements
# Ryo Takatori
# 10/16/2020

# Modules
import sys
sys.path.append("orbit_tools")
import numpy as np
import coord_trans

# State Vector to Classical Orbital Elements
def rv2coe(mu, r, v, deg):
    # Inertial position and velocity
    r = np.array(r)
    v = np.array(v)
    norm_r = np.linalg.norm(r) # position magnitude
    norm_v = np.linalg.norm(v) # velocity magnitude
 
    # Semi-major axis (a)
    epsilon = ((norm_v**2)/2) - (mu/norm_r) # specific mechanical energy [MJ/kg]
    a = -mu/(2*epsilon)

    # Eccentricity (ecc)
    e_vec = (1/mu)*((((norm_v**2)-(mu/norm_r))*r)-(np.dot(r,v)*v)) # eccentricity vector
    ecc = np.linalg.norm(e_vec)

    # Inclination (inc)
    h_vec = np.cross(r,v)
    h = np.linalg.norm(h_vec)
    inc = np.arccos(h_vec[2]/h)

    # Right Ascension of Ascending Node (raan)
    k_hat = np.array([0, 0, 1])
    i_hat = np.array([1, 0, 0])
    n = np.cross(k_hat,h_vec)
    norm_n = np.linalg.norm(n)
    raan = np.arccos(np.dot(i_hat,n)/norm_n)
    # Quadrant ambiguity check
    if n[1] < 0:
        raan = 2*np.pi - raan
    
    # Argument of perigee (w)
    w = np.arccos(np.dot(n,e_vec)/(norm_n*ecc))
    # Quadrant ambiguity check
    if e_vec[2] < 0:
        w = 2*np.pi - w

    # True Anomaly (theta)
    theta = np.arccos(np.dot(e_vec,r)/(ecc*norm_r))
    # Quadrant ambiguity check
    if np.dot(r,v) < 0:
        theta = 2*np.pi - theta
    
    # Classical orbital elements
    coe = [a, ecc, h, inc, raan, w, theta]

    # Convert to degrees if needed
    if deg:
        coe[3:] = np.rad2deg(coe[3:])
    return coe

# Classical Orbital Elements to State Vector
def coe2rv(mu, coe, deg):
    if deg:
        coe[3:] = np.deg2rad(coe[3:])
    [a, ecc, h, inc, raan, w, theta] = coe
    # State vectors in perifocal
    r_peri = ((h**2)/mu)*(1/(1+ecc*np.cos(theta)))*np.array([np.cos(theta), np.sin(theta), 0])
    v_peri = (mu/h)*np.array([-np.sin(theta), ecc+np.cos(theta), 0])
    # Convert
    [r,v] = coord_trans.peri2eci(r_peri,v_peri,w,raan,inc)
    return r,v
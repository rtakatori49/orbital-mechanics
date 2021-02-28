# Equation of Motion
# Ryo Takatori
# 10/16/2020

# Modules
import numpy as np
import celestial_body_ephemeris as cbe
import perturbations
import planet_data as pd
import time_trans as tt

# Two body
def two_body(t, y, mu):
    rx, ry, rz, vx, vy, vz = y # velocity
    r = np.array([rx,ry,rz]) # position
    norm_r = np.linalg.norm(r) # position magnitude
    ax, ay, az = -r*mu/norm_r**3 # acceleration
    return [vx, vy, vz, ax, ay, az]

def cowell(t, y, drag, srp, n_body, j2_6, mu, A, m, c_d, c_r, epoch, planet):
    rx, ry, rz, vx, vy, vz = y # velocity
    r = np.array([rx,ry,rz]) # position
    norm_r = np.linalg.norm(r) # position magnitude
    v = np.array([vx,vy,vz])
    a_zero = np.array([0,0,0])
    if drag:
        a_drag = perturbations.drag_exp(r, v, A, m, c_d)
    else:
        a_drag = a_zero
    if srp:
        a_srp = perturbations.srp(t, epoch, r, A, m, c_r)
    else:
        a_srp = a_zero
    if n_body:
        if planet == "Moon":
            jd = tt.jd(epoch)
            r_body = cbe.moon(jd)
            mu_body = pd.Moon.mu
        a_n_body = perturbations.n_body(t, r, r_body, mu_body)
    else:
        a_n_body = a_zero
    if j2_6:
        a_j2_6 = perturbations.j2_6(r)
    else:
        a_j2_6 = a_zero
    ax, ay, az = -r*mu/norm_r**3 + a_drag + a_srp + a_n_body + a_j2_6 # acceleration
    return [vx, vy, vz, ax, ay, az]
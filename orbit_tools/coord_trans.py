# Coordinate Tranformation
# Ryo Takatori
# 10/16/2020

# Modules
import sys
sys.path.append("orbit_tools")
import numpy as np
import dcm

# Perifocal to Earth Centered Inertial (ECI)
def peri2eci(r_peri,v_peri,w,raan,inc):
    # Rotation 3 through argument of perigee
    cz_w = dcm.rotz(w)
    cz_w[1,0] = -cz_w[1,0]
    cz_w[0,1] = -cz_w[0,1]

    # Rotation 1 through inclination
    cx_inc = dcm.rotx(inc)
    cx_inc[1,2] = -cx_inc[1,2]
    cx_inc[2,1] = -cx_inc[2,1]

    # Rotation 3 through RAAN
    cz_raan = dcm.rotz(raan)
    cz_raan[1,0] = -cz_raan[1,0]
    cz_raan[0,1] = -cz_raan[0,1]

    # DCM
    Q = np.matmul(cz_w,np.matmul(cx_inc,cz_raan))
    Q = np.transpose(Q)

    # States
    r = np.matmul(Q,np.transpose(r_peri))
    v = np.matmul(Q,np.transpose(v_peri))
    return r,v

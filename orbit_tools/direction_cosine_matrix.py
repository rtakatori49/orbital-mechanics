# Direction Cosine Matrix
# Ryo Takatori
# 10/16/20202

# Modules
import numpy as np

# Rotation about x
def rot_x(theta):
    cx = np.array([[1, 0, 0],\
        [0, np.cos(theta), -np.sin(theta)],\
        [0, np.sin(theta), np.cos(theta)]])
    return cx

# Rotation about y
def rot_y(theta):
    cy = np.array([[np.cos(theta), 0, np.sin(theta)],\
        [0, 1, 0],\
        [-np.sin(theta), 0, np.cos(theta)]])
    return cy

# Rotation about z    
def rot_z(theta):
    cz = np.array([[np.cos(theta), -np.sin(theta), 0],\
        [np.sin(theta), np.cos(theta), 0],\
        [0, 0, 1]])
    return cz

# Clohessy–Wiltshire Matrix
def clohessy_wiltshire(n, t):
    nt = n*t
    phi_rr = np.array([[4-3*np.cos(nt), 0, 0],\
        [6*(np.sin(nt)-nt), 1, 0],\
        [0, 0, np.cos(nt)]])
    phi_rv = np.array([[(1/n)*np.sin(nt), (2/n)*(1-np.cos(nt)), 0],\
        [(2/n)*(np.cos(nt)-1), (1/n)*(4*np.sin(nt)-3*(nt)), 0],\
        [0, 0, (1/n)*np.sin(nt)]])
    phi_vr = np.array([[3*n*np.sin(nt), 0, 0],\
        [6*n*(np.cos(nt)-1), 0, 0],\
        [0, 0, -n*np.sin(nt)]])
    phi_vv = np.array([[np.cos(nt), 2*np.sin(nt), 0],\
        [-2*np.sin(nt), 4*np.cos(nt)-3, 0],\
        [0, 0, np.cos(nt)]])
    return phi_rr, phi_rv, phi_vr, phi_vv
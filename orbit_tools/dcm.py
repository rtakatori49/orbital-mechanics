# Direction Cosine Matrix
# Ryo Takatori
# 10/16/20202

# Modules
import numpy as np

# Rotation about x
def rotx(theta):
    cx = [[1, 0, 0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]]
    cx = np.array(cx)
    return cx

# Rotation about y
def roty(theta):
    cy = [[np.cos(theta), 0, np.sin(theta)], [0, 1, 0], [-np.sin(theta), 0, np.cos(theta)]]
    cy = np.array(cy)
    return cy

# Rotation about z    
def rotz(theta):
    cz = [[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]]
    cz = np.array(cz)
    return cz

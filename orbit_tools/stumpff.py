# Stumpff Functions
import numpy as np

# Stump funtion C by Curtis
def c(z):
    # This function evaluates the Stumpff function C(z) according
    # to Equation 3.53.
    # z - input argument
    # c - value of C(z)
    if z > 0:
        c = (1 - np.cos(np.sqrt(z)))/z
    elif z < 0:
        c = (np.cosh(np.sqrt(-z)) - 1)/(-z)
    else:
        c = 1/2
    return c

# Stump function S by Curtis
def s(z):
    # This function evaluates the Stumpff function S(z) according
    # to Equation 3.52.
    # z - input argument
    # s - value of S(z)
    if z > 0:
        s = (np.sqrt(z) - np.sin(np.sqrt(z)))/(np.sqrt(z))**3
    elif z < 0:
        s = (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/(np.sqrt(-z))**3
    else:
        s = 1/6
    return s
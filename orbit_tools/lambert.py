
import numpy as np
from . import stumpff
# Lamberts Solution

# Universal Variable
def universal_variable(r0, rf, dt, mu, tm=1):
    r0 = np.array(r0)
    rf = np.array(rf)
    norm_r0 = np.linalg.norm(r0)
    norm_rf = np.linalg.norm(rf)
    cross_r = np.cross(r0, rf)
    dtheta = np.arccos(np.dot(r0,rf)/(norm_r0*norm_rf))
    if tm == 1:
        if cross_r[2] < 0:
            dtheta = 2*np.pi - dtheta
    elif tm == -1:
        if cross_r[2] >= 0:
            dtheta = 2*np.pi - dtheta
    A = np.sin(dtheta)*np.sqrt(norm_r0*norm_rf/(1-np.cos(dtheta)))

    # Functions
    def C(z):
        return stumpff.c(z)
    def S(z):
        return stumpff.s(z)
    def y(z):
        dum = norm_r0 + norm_rf + A*(z*S(z)-1)/np.sqrt(C(z))
        return dum
    def F(z, t):
        if y(z) < 0:
            dum = float("nan")
        else:
            dum = ((y(z)/C(z))**1.5)*S(z) + A*np.sqrt(y(z)) - np.sqrt(mu)*t
        return dum
    def dFdz(z):
        if z == 0:
            dum = np.sqrt(2)/40*y(0)**1.5 + A/8*(np.sqrt(y(0)) + A*np.sqrt(1/2/y(0)))
        else:
            dum = ((y(z)/C(z))**1.5)*(1/2/z*(C(z) - 3*S(z)/2/C(z)) \
                + 3*S(z)**2/4/C(z)) + A/8*(3*S(z)/C(z)*np.sqrt(y(z)) \
                + A*np.sqrt(C(z)/y(z)))
        return dum
    z = -100
    while F(z, dt) < 0 or np.isnan(F(z, dt)):
        z += 0.1
    tol = 10e-8
    ratio = 1
    n = 0
    nmax = 5000
    while np.abs(ratio) > tol and n < nmax:
        n += 1
        ratio = F(z, dt)/dFdz(z)
        z -= ratio
    f = 1 - y(z)/norm_r0
    g = A*np.sqrt(y(z)/mu)
    gdot = 1 - y(z)/norm_rf
    v0 = (1/g)*(rf-f*r0)
    vf = (1/g)*(gdot*rf-r0)
    return v0, vf

# Gauss
def gauss(r0, rf, dt, tm, mu):
    r0 = np.array(r0)
    rf = np.array(rf)
    norm_r0 = np.linalg.norm(r0)
    norm_rf = np.linalg.norm(rf)
    cosdtheta = np.dot(r0,rf)/(norm_r0*norm_rf)
    dtheta = np.arccos(cosdtheta)
    cosdtheta_half = np.sqrt((1+cosdtheta)/2)
    sindtheta = tm*np.sqrt(1-cosdtheta**2)
    # l and m calculation
    L = ((norm_rf+norm_r0)/(4*np.sqrt(norm_rf*norm_r0)*cosdtheta_half)-0.5)
    m = (mu*dt**2)/(2*np.sqrt(norm_rf*norm_r0)*cosdtheta_half)**3
    # Iterate for y and find x1
    # Initialize
    y = [1]
    tol = 1e-6
    diff = 1000
    j = 0
    while diff > tol:
        if diff > 10000:
            break
        x1 = (m/(y[j]**2)) - L
        # Series expansion
        x2 = (4/3)*(1+((6*x1)/5)+((6*8*x1**2)/(5*7))+((6*8*10*x1**3)/(5*7*9)))
        y.append(1+x2*(L+x1))
        # Check differences
        diff = np.abs(y[j]-y[j+1])
        j += 1
    # Calculate parameter
    p = (norm_r0*norm_rf*(1-cosdtheta))/(norm_r0+norm_rf-2*np.sqrt(norm_r0*norm_rf)*cosdtheta_half*(1-(2*x1)))

    # Lagrange coefficients
    f = 1 - (norm_rf/p)*(1-cosdtheta)
    g = (norm_r0*norm_rf*sindtheta)/(np.sqrt(mu*p))
    fdot = np.sqrt(1/p)*np.tan(dtheta/2)*(((1-cosdtheta)/p)-(1/norm_rf)-(1/norm_r0))
    gdot = 1 - (norm_r0/p)*(1-cosdtheta)
    # Velocity vectors [km/s]
    v0 = (rf-f*r0)/g
    vf = (gdot*rf-r0)/g
    return v0, vf

# Izzo Gooding
def izzo_gooding():
    pass

# Minimum Energy
def minimum_energy(r0, rf, tm, mu):
    r0 = np.array(r0)
    rf = np.array(rf)
    norm_r0 = np.linalg.norm(r0)
    norm_rf = np.linalg.norm(rf)
    # Transfer angle [rad]
    dtheta = np.arccos(np.dot(r0,rf)/(norm_r0*norm_rf))
    dtheta = np.arcsin(tm*np.sqrt(1-np.cos(dtheta)**2))
    # Chord length
    c = np.sqrt(norm_r0**2+norm_rf**2-2*norm_r0*norm_rf*np.cos(dtheta))
    # Semi-perimeter
    s = (norm_r0+norm_rf+c)/2
    # Minimum semi-major axis
    amin = s/2
    # Minimum parameter
    pmin = (norm_r0*norm_rf*(1-np.cos(dtheta)))/c
    # Minimum eccentricity
    emin = np.sqrt(1-((2*pmin)/s))
    # Alpha and beta
    alphae = np.pi
    betae = 2*np.arcsin(np.sqrt((s-c)/s))
    # Minimum energy transfer time
    dtmin = np.sqrt(amin**3/mu)*(alphae-tm*(betae-np.sin(betae)))
    dtp = (1/3)*np.sqrt(2/mu)*((s**(3/2))-(s-c)**(3/2))
    # Velocity vector [km/s]
    v0 = (np.sqrt(mu*pmin)/(norm_r0*norm_rf*np.sin(dtheta)))*(rf-(1-(norm_rf/pmin)*(1-np.cos(dtheta)))*r0)
    return amin, emin, dtmin, dtp, v0
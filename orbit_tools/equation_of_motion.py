# Equation of Motion
# Ryo Takatori
# 10/16/2020

# Modules
import numpy as np
from . import classical_orbital_elements as coe
from . import coordinate_transformation as c_t
from . import direction_cosine_matrix as dcm
from . import perturbations as pert
from . import stumpff

# Two body
def two_body(t, y, mu):
    rx, ry, rz, vx, vy, vz = y # velocity
    r = np.array([rx,ry,rz]) # position
    norm_r = np.linalg.norm(r) # position magnitude
    ax, ay, az = -r*mu/norm_r**3 # acceleration
    return [vx, vy, vz, ax, ay, az]

# Altitude Checker:
# Cowell
def alt_checker(t, y, mu, perts):
    return np.linalg.norm(y[0:3]) - (perts.radius+perts.alt)
alt_checker.terminal = True
alt_checker.direction = -1
# VoP
def alt_checker_vop(t, y, mu, perts):
    r, _ = coe.coe_to_state(mu, y, deg=False)
    return np.linalg.norm(r) - (perts.radius+perts.alt)
alt_checker_vop.terminal = True
alt_checker_vop.direction = -1

# Cowell"s Method
def cowell(t, y, mu, perts):
    rx, ry, rz, vx, vy, vz = y # velocity
    r = np.array([rx,ry,rz]) # position
    norm_r = np.linalg.norm(r) # position magnitude
    v = np.array([vx,vy,vz])
    a_pert = pert.acceleration(t, r, v, perts) # perturbations
    ax, ay, az = -r*mu/norm_r**3 + a_pert # acceleration
    return [vx, vy, vz, ax, ay, az]

# Variation of Parameters
def vop(t, y, mu, perts):
    # Get vectors
    [ecc, h, inc, _, w, theta] = y
    r, v = coe.coe_to_state(mu, y, deg=False) # vectors
    norm_r = np.linalg.norm(r) # position magnitude
    rv_cross = np.cross(r,v)

    # Perturbing acceleration
    n_hat = rv_cross/np.linalg.norm(rv_cross)
    r_hat = r/norm_r
    nr_cross = np.cross(n_hat, r)
    t_hat = nr_cross/np.linalg.norm(nr_cross)
    a_pert = pert.acceleration(t, r, v, perts) # perturbations
    R = np.dot(a_pert, r_hat)
    T = np.dot(a_pert, t_hat)
    N = np.dot(a_pert, n_hat)
    
    # d(element)/dt calculation
    # Angular momentum
    dhdt = norm_r*T
    
    # Eccentricity
    dedt = (h/mu)*np.sin(theta)*R \
            + (1/(mu*h))*((h**2+(mu*norm_r))*np.cos(theta)+mu*ecc*norm_r)*T
    
    # True anomaly
    dthetadt_two_body = h/norm_r**2
    dthetadt_pert = (1/(ecc*h))*((((h**2)*R)/mu)*np.cos(theta) \
        - (((h**2)/mu)+norm_r)*T*np.sin(theta))
    dthetadt_all = dthetadt_two_body + dthetadt_pert

    # Inclination
    u = w + theta
    dincdt = (norm_r/h)*N*np.cos(u)

    # Right ascention of ascending node
    dRAANdt = ((norm_r*np.sin(u))/(h*np.sin(inc)))*N

    # Argument of periapse
    dwt = ((-norm_r*np.sin(u))/(h*np.tan(inc)))*N - dthetadt_pert

    # Return states
    return [dedt, dhdt, dincdt, dRAANdt, dwt, dthetadt_all]

# Universal Variable by Curtis
def kepler_u(d_t, norm_r0, r0_v0, alpha, mu):
    # This function uses Newton"s method to solve the universal
    # Kepler equation for the universal anomaly.
    # mu - gravitational parameter (km^3/s^2)
    # x - the universal anomaly (km^0.5)
    # dt - time since x = 0 (s)
    # ro - radial position (km) when x = 0
    # vro - radial velocity (km/s) when x = 0
    # a - reciprocal of the semimajor axis (1/km)
    # z - auxiliary variable (z = a*x^2)
    # C - value of Stumpff function C(z)
    # S - value of Stumpff function S(z)
    # n - number of iterations for convergence
    # nMax - maximum allowable number of iterations
    error = 1.0e-8
    n_max = 1000
    x = np.sqrt(mu)*np.abs(alpha)*d_t
    n = 0
    ratio = 1
    while np.abs(ratio) > error and n <= n_max:
        n = n + 1
        c = stumpff.c(alpha*x**2)
        s = stumpff.s(alpha*x**2)
        F = (norm_r0*r0_v0*x**2*c)/np.sqrt(mu) \
            + (1-alpha*norm_r0)*x**3*s + norm_r0*x - np.sqrt(mu)*d_t
        dFdx = (norm_r0*r0_v0*x*(1-alpha*x**2*s))/np.sqrt(mu) \
            + (1-alpha*norm_r0)*x**2*c + norm_r0
        ratio = F/dFdx
        x = x - ratio
    if n > n_max:
        print(f"")
    return x

# Lagrange F and G coefficients by Curtis
def f_and_g(x, t, norm_r0, alpha, mu):
    # This function calculates the Lagrange f and g coefficients.
    # mu - the gravitational parameter (km^3/s^2)
    # a - reciprocal of the semimajor axis (1/km)
    # ro - the radial position at time to (km)
    # t - the time elapsed since ro (s)
    # x - the universal anomaly after time t (km^0.5)
    # f - the Lagrange f coefficient (dimensionless)
    # g - the Lagrange g coefficient (s)
    z = alpha*x**2
    f = 1 - (x**2/norm_r0)*stumpff.c(z)
    g = t - (x**3*stumpff.s(z))/np.sqrt(mu)
    return f, g

# Largrange F and G derivative coefficients by Curtis
def fdot_and_gdot(x, norm_r, norm_r0, alpha, mu):
    # This function calculates the time derivatives of the
    # Lagrange f and g coefficients.
    # mu - the gravitational parameter (km^3/s^2)
    # a - reciprocal of the semimajor axis (1/km)
    # ro - the radial position at time to (km)
    # t - the time elapsed since initial state vector (s)
    # r - the radial position after time t (km)
    # x - the universal anomaly after time t (km^0.5)
    # fdot - time derivative of the Lagrange f coefficient (1/s)
    # gdot - time derivative of the Lagrange g coefficient (dimensionless)
    z = alpha*x**2
    fdot = (np.sqrt(mu)*(z*stumpff.s(z)-1)*x)/(norm_r*norm_r0)
    gdot = 1 - (x**2*stumpff.c(z))/norm_r
    return fdot, gdot

# Kepler"s Equation based on Curtis
def kepler(r0, v0, t, mu):
    # This function computes the state vector (r,v) from the
    # initial state vector (r0,v0) and the elapsed time.
    # mu - gravitational parameter (km^3/s^2)
    # r0 - initial position vector (km)
    # v0 - initial velocity vector (km/s)
    # t - elapsed time (s)
    # r - final position vector (km)
    # v - final velocity vector (km/s)
    
    # Magnitudes of r0 and v0
    r0 = np.array(r0)
    v0 = np.array(v0)
    norm_r0 = np.linalg.norm(r0)
    norm_v0 = np.linalg.norm(v0)
    r0_v0 = np.dot(r0, v0)/norm_r0 # initial radial velocity
    alpha = (2/norm_r0) - ((norm_v0**2)/mu) # reciprocal of the semi-major axis
    x = kepler_u(t, norm_r0, r0_v0, alpha, mu) # universal anomaly
    f, g = f_and_g(x, t, norm_r0, alpha, mu) # f and g functions
    r = f*r0 + g*v0 # final position vector
    norm_r = np.linalg.norm(r) # magnitude of r
    fdot, gdot = fdot_and_gdot(x, norm_r, norm_r0, alpha, mu) # derivatives of da and g
    v = fdot*r0 + gdot*v0 # final velocity vector
    return r, v

# Encke"s Method
def encke(t_span, d_t, r, v, mu, tol, perts):
    def parse_return(r, v, t):
        y = np.append(r, v, axis=1)
        y = np.transpose(y)
        time_list = [num/(24*60*60) for num in [*range(0, t, d_t)] + [t + d_t]]
        return y, time_list

    # Set osculating vectors to current vectors
    r_osc = r # inital position as osculating orbit position [km]
    v_osc = v # inital velocity as osculating orbit velocity [km/s]
    # Append to array
    r = np.array([r])
    v = np.array([v])
    # Set relative vectors to zero
    zero_array = np.array([0,0,0])
    d_r = zero_array # relative position of perturbing orbit position [km]
    d_v = zero_array # relative velocity of perutbing orbit velocity [km/s]
    for x, t in enumerate(range(0, t_span[1], d_t)):
        r_osc, v_osc = kepler(r_osc, v_osc, d_t, mu) # propgate osculating orbit using Universal Variable 
        # Calculate binomial expansion values (q & F)
        q = np.dot(d_r,(2*r[x,:]-d_r))/np.linalg.norm(r[x,:])**2
        F = (((q**2)-(3*q)+3)/(1+((1-q)**(3/2))))*q
        norm_r_osc = np.linalg.norm(r_osc) # magnitude of osculating orbit position [km]
        a_pert = pert.acceleration(t, r[x,:], v[x,:], perts) # perturbations
        d_a = a_pert - (mu/norm_r_osc**3)*(d_r-F*r[x,:]) # relative acceleration of perturbing orbit acceleration [km/s^2]
        d_v = d_a*d_t + d_v # relative velocity of perutbing orbit velocity [km/s]
        d_r = 0.5*d_a*(d_t**2) + d_v*d_t + d_r # relative position of perturbing orbit position [km]
        # Check for altitude
        if np.linalg.norm(r[x,:]) < perts.radius + perts.alt:
            print("Altitude too low")
            y, time_list = parse_return(r, v, t)
            return r, v, y, time_list
        # Calculate actual vectors
        r = np.append(r, [r_osc+d_r], axis=0) # actual position vector [km]
        v = np.append(v, [v_osc+d_v], axis=0) # actual velocity vector [km/s]
        # Check for rectification Rectify if true
        if np.linalg.norm(d_r)/np.linalg.norm(r[x,:])>tol:
            r_osc = r[x+1,:] # actual position as osculating orbit position [km]
            v_osc = v[x+1,:] # actual velocity as osculating orbit velocity [km/s]
            d_r = zero_array # relative position of perturbing orbit position [km]
            d_v = zero_array # relative velocity of perutbing orbit velocity [km/s]
    y, time_list = parse_return(r, v, t)
    return r, v, y, time_list

# Kepler's Equation
def kepler_e(ecc, me):
    if me < np.pi:
        E = me + ecc/2
    else:
        E = me - ecc/2
    ratio = 2
    error = 10^-8
    i = 0
    while error < abs(ratio):
        ratio = (me-E+ecc*np.sin(E))/(1-ecc*np.cos(E))
        E = E + ratio
        i += 1
        if i == 10000:
            break
    return E

# Circular Restricted Three Body Problem
def cr3bp(t, y, mustar):
    # Position & Velocity
    rx, ry, vx, vy = y
    r_1 = np.sqrt(((rx-mustar)**2)+(ry**2))
    r_2 = np.sqrt(((rx+1-mustar)**2)+(ry**2))

    # Acceleration
    ax = rx + (2*vy) - (((1-mustar)*(rx-mustar))/r_1**3) - ((mustar*(rx+1-mustar))/r_2**3)
    ay = ry - (2*vx)- (((1-mustar)*ry)/r_1**3) - ((mustar*ry)/r_2**3)
    return [vx, vy, ax, ay]

# Five Term Acceleration
def fta(r_A, r_B, v_A, v_B, a_A, a_B):
    # Coordinate transformation
    _, q = c_t.eci_to_lvlh(r_A, v_A) # ECI to LVLH DCM

    # Relative position
    r_B_A = r_B - r_A # relative position [km]
    r_B_A_LVLH = np.dot(q,r_B_A) # position in LVLH [km]

    # Relative velocity
    norm_r_A = np.linalg.norm(r_A) # position magnitude [km]
    h_A = np.cross(r_A, v_A) # angular momentum [km^2/s]
    omega = h_A/norm_r_A**2 # relative angular velocity [rad/s]
    v_B_A = v_B - v_A - np.cross(omega, r_B_A) # relative velocity [km/s]
    v_B_A_LVLH = np.dot(q, v_B_A) # relative velocity in LVLH [km/s]

    # Relative acceleration
    omega_dot = ((-2*np.dot(v_A, r_A))*omega)/norm_r_A**2 # Relative angular acceleration [rad/s^2]
    a_B_A = a_B - a_A - np.cross(omega_dot, r_B_A) - np.cross(omega,np.cross(omega,r_B_A)) - 2*np.cross(omega,v_B_A) # Relative acceleration [km/s^2]
    a_B_A_LVLH = np.dot(q, a_B_A) # Relative acceleration in LVLH [km/s^2]
    return r_B_A_LVLH, v_B_A_LVLH, a_B_A_LVLH

# Rendezvous maneuvers

# Two Impulse Maneuver
def two_impluse(r_rel, v_rel, n, t):
    # Matrices
    phi_rr, phi_rv, phi_vr, phi_vv = dcm.clohessy_wiltshire(n, t)
    # Velocity to get on [km/s]
    v_plus_2 = np.dot(np.dot(np.linalg.inv(phi_rv), (-phi_rr)), r_rel)
    delta_v_1 = v_plus_2 - v_rel
    # Velocity to get off [km/s]
    v_f_minus_2 = np.dot(phi_vr, r_rel) + np.dot(phi_vv, v_plus_2)
    delta_v_2 = -v_f_minus_2
    # Delta-V [km/s]
    delta_v = np.linalg.norm(delta_v_1) + np.linalg.norm(delta_v_2)
    return delta_v, delta_v_1, delta_v_2
    
# Linearized Equations of Motion
def leom(t, y, mu, r_bar, v_bar, v_c, t0):
    # Unpack
    rx, ry, rz, vx, vy, vz, del_x, del_y, del_z, del_x_dot, del_y_dot, del_z_dot = y # velocity
    ## Chief
    r = np.array([rx,ry,rz]) # position
    norm_r = np.linalg.norm(r) # position magnitude
    v = np.array([vx,vy,vz]) # position
    h = np.linalg.norm(np.cross(r, v))
    ax, ay, az = -r*mu/norm_r**3 # acceleration

    # Deputy
    del_x_2dot = (((2*mu)/norm_r**3)+(h**2/norm_r**4))*del_x - ((2*(np.dot(v,r))*h)/norm_r**4)*del_y+(2*h/norm_r**2)*del_y_dot
    if r_bar:
        del_x_2dot -= (3*mu*del_x)/norm_r**3
    if v_bar:
        del_x_2dot -= 2*np.sqrt(mu/norm_r**3)*v_c
    del_y_2dot = -(((mu)/norm_r**3)-(h**2/norm_r**4))*del_y + ((2*(np.dot(v,r))*h)/norm_r**4)*del_x-(2*h/norm_r**2)*del_x_dot
    if v_bar and t == t0:
        del_y_2dot += v_c*t
    del_z_2dot = -(mu/norm_r**3)*del_z
    return [vx, vy, vz, ax, ay, az, del_x_dot, del_y_dot, del_z_dot, del_x_2dot, del_y_2dot, del_z_2dot]

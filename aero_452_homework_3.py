# Homework 3
# AERO 452
# Ryo Takatori

import numpy as np
from orbit_tools import classical_orbital_elements as coe
from orbit_tools import equation_of_motion as eom
from orbit_tools import perturbations
from orbit_tools import planet_data as p_d
from orbit_tools import plotter
from scipy.integrate import solve_ivp
import time

earth = p_d.Earth()
mue = earth.mu
re = earth.radius

def problem_1():
    # Problem 1
    print("Problem 1:")

    t_total = 120*24*60*60 # Time [days=>s]
    t_span = [0, t_total]

    # Spacecarft
    d = 1 # Diameter [m]    
    mass = 100 # Mass [kg]
    area = np.pi*(d/2)**2 # Ram area [m^2]
    c_d = 2.2 # Coefficient of drag

    # Initial Orbit
    zp = 215 # Perigee altitude [km]
    za = 939 # Apogee altitude [km]
    raan = 340 # RAAN [deg]
    inc= 65.2 # Inclination [deg]
    w = 58 # Arguement of perigee [deg]
    theta = 332 # True Anomaly [deg]
    rp = re + zp # Radius of perigee [km]
    ra = re + za # Radius of apogee [km]
    a = (rp+ra)/2 # Semi-major axis [km]
    ecc = (ra-rp)/(ra+rp) # Eccentricity
    h = np.sqrt(a*mue*(1-ecc**2)) # Angular momentum [km^2/s]
    coe_ = [a, ecc, h, inc, raan, w, theta]
    perts = perturbations.Perturbation()
    perts.set_drag(area, mass, c_d)

    # Cowell
    r0, v0 = coe.coe_to_state(mue, coe_[1:7])
    y0 = np.concatenate((r0, v0))
    tic = time.time()
    sol = solve_ivp(eom.cowell, t_span, y0, method="RK45", events=eom.alt_checker, args=(mue, perts,), rtol=1e-8, atol=1e-8)
    toc = time.time()
    print(f"(Cowell) : Elapsed Time {str(toc-tic)} [s]")
    plotter.coe_plot(sol.y[0:6], sol.t/(24*60*60), "vector", mue, re, "Problem 1 Cowell", deg=True, diff=False)

    # VoP
    y0_vop = coe_[1:3] + list(np.deg2rad(coe_[3:7]))
    tic = time.time()
    sol = solve_ivp(eom.vop, t_span, y0_vop, method="RK45", events=eom.alt_checker_vop, args=(mue, perts,), rtol=1e-8, atol=1e-8)
    toc = time.time()
    print(f"(VoP) : Elapsed Time {str(toc-tic)} [s]")
    plotter.coe_plot(sol.y[0:6], sol.t/(24*60*60), "coe", mue, re, "Problem 1 VoP", deg=False, diff=False)

    # Encke
    d_t = 600
    tic = time.time()
    _, _, y, t = eom.encke(t_span, d_t, r0, v0, mue, 1e-8, perts)
    toc = time.time()
    print(f"(Encke) : Elapsed Time {str(toc-tic)} [s]")
    plotter.coe_plot(y, t, "vector", mue, re, "Problem 1 Encke", deg=True, diff=False)

    # Two-Body
    tic = time.time()
    sol = solve_ivp(eom.two_body, t_span, y0, method="RK45", args=(mue,), rtol=1e-8, atol=1e-8)
    toc = time.time()
    print(f"(Two-Body) : Elapsed Time {str(toc-tic)} [s]")
    plotter.coe_plot(sol.y[0:6], sol.t/(24*60*60), "vector", mue, re, "Problem 1 Two-Body", deg=True, diff=False)

    print("Unintuitively, Encke's method takes the longest. This mostly "\
        "has to do with the fact that the time step at which I am stepping is "\
        " far smaller than the steps that the ode45 is stepping at. However, if I "\
        "try to step at a similar interval to the ode45, Encke's gives a different "\
        "plot. I assumed that this was because the perturbational effects of "\
        "atmospheric drag is significant at this altitude. Furthermore, since "\
        "the perigee is so close to the altitude at which we consider the orbit "\
        "to be no longer space, if Encke and Universal Variable cuts the corner "\
        "past the 100 [km] point it will prematurely return the code. The density "\
        "is also rapidly pulling the spacecraft down since its simulating a "\
        "density that is not accurate which further accelerates the spacecraft "\
        "down towards Earth. The other methods were also interesting because "\
        "the Cowell""s method actually took longer than the VoP. This was "\
        "unexpected as I had thought that Cowell""s would do better with less terms. "\
        "The two-body problem is as expected, there are no changes over time. "\
        "There are barely any changes in the COE as well.")

def problem_2():
    ## Problem 2
    print("Problem 2:")

    t_total = 48*60*60 # Time [hr=>s]
    t_span = [0, t_total]

    zp = 300 # Perigee altitude [km]
    za = 3092 # Apogee altitude [km]
    raan = 45 # RAAN [deg]
    inc = 28 # Inclination [deg]
    w = 30 # Arguement of perigee [deg]
    theta = 40 # True Anomaly [deg]
    rp = re + zp # Radius of perigee [km]
    ra = re + za # Radius of apogee [km]
    a = (rp+ra)/2 # Semi-major axis [km]
    ecc = (ra-rp)/(ra+rp) # Eccentricity
    h = np.sqrt(a*mue*(1-ecc**2)) # Angular momentum [km^2/s]
    coe_ = [a, ecc, h, inc, raan, w, theta]
    r0, v0 = coe.coe_to_state(mue, coe_[1:7])
    y0 = np.concatenate((r0, v0))

    # J2
    perts = perturbations.Perturbation()
    perts.set_j2_6(2)
    tic = time.time()
    sol = solve_ivp(eom.cowell, t_span, y0, method="RK45", events=eom.alt_checker, args=(mue, perts,), rtol=1e-8, atol=1e-8)
    toc = time.time()
    print(f"(Cowell J2) : Elapsed Time {str(toc-tic)} [s]")
    plotter.coe_plot(sol.y[0:6], sol.t/(24*60*60), "vector", mue, re, "Problem 2 Cowell J2", deg=True, diff=False)

    # J3
    perts = perturbations.Perturbation()
    perts.set_j2_6(3)
    tic = time.time()
    sol = solve_ivp(eom.cowell, t_span, y0, method="RK45", events=eom.alt_checker, args=(mue, perts,), rtol=1e-8, atol=1e-8)
    toc = time.time()
    print(f"(Cowell J3) : Elapsed Time {str(toc-tic)} [s]")
    plotter.coe_plot(sol.y[0:6], sol.t/(24*60*60), "vector", mue, re, "Problem 2 Cowell J3", deg=True, diff=False)

    print("This is the expected plots as the J2 should have a significant "\
        "effect over the J3. We can see in the plots that they are nearly "\
        "identical. We went over in class that the effect of J2 should be three "\
        "orders of magnitude higher, which means that the plots should barely be "\
        "changing when added together.")

def problem_3():
    ## Problem 3
    print("Problem 3:")


def main():
    problem_1()
    problem_2()
    problem_3()

if __name__=="__main__":
    main()
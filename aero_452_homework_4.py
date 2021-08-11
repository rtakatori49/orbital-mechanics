# Homework 4
# AERO 452
# Ryo Takatori

import datetime
import matplotlib.pyplot as plt
import numpy as np
from orbit_tools import equation_of_motion as eom
from orbit_tools import perturbations
from orbit_tools import planet_data as p_d
from orbit_tools import plotter
from orbit_tools import time_trans as t_t
from orbit_tools import utility as util
from scipy.integrate import solve_ivp
import time

earth = p_d.Earth()
mue = earth.mu
re = earth.radius
moon = p_d.Moon()

def problem_1():
    # Problem 1
    print("Problem 1:")
    start_date = datetime.datetime(1997, 8, 10, 0, 0, 0)
    end_date = datetime.datetime(1997, 8, 24, 0, 0, 0)
    t_duration = (end_date - start_date).total_seconds()
    r0 = np.array([-26175.1034, 12757.0706, 14626.6556]) # Position vector [km]
    v0 = np.array([2.376641, 0.139677, 2.078097]) # Velocity vector [km/s]
    t_span = [0, t_duration]
    y0 = np.concatenate((r0, v0))
    sol = solve_ivp(eom.two_body, t_span, y0, method="RK45", args=(mue,), rtol=1e-8, atol=1e-8)
    jd = t_t.jd(start_date) # Julian date
    t_prev = t_total = 0
    nu_total = []
    for x, y, z, t, in zip(sol.y[0], sol.y[1], sol.y[2], sol.t):
        r = np.array([x, y, z])
        nu = perturbations.shadow(jd, t, r) # check shadow
        nu_total.append(nu)
        if nu == 1:
            t_diff = t - t_prev
            t_total += t_diff
        t_prev = t
    eclipse_time = (t_total-t_duration)/(24*60*60)
    percent_lit = (t_total/t_duration)*100

    # Display
    print(f"Time spacecraft is in Eclipse: {eclipse_time} [day]")
    print(f"Percentage of time spacecraft is lit: {percent_lit} [%]")
    print(f"Percent increase of solar radiation pressure with assumption of continuous Sunlight: {100-percent_lit} [%]")

    plt.plot([t/(24*60*60) for t in sol.t], nu_total)
    plt.title("Shadow Over Time")
    plt.xlabel("Time [day]")
    plt.ylabel("Light")
    plt.ylim([-0.5, 1.5])
    plt.show()

    print("Over the two week period, the spacecraft spent most of its "\
        "time in the Sun. There were less than 5 % time in eclipse. This "\
        "makes sense because the orbit was relatively big, so it was spending "\
        "more time exposed to the sun. With this small variation, ignoring it "\
        "would show some difference, but assuming SRP to be on all the time "\
        "would not be a terrible assumption.")

def problem_2():
    ## Problem 2
    print("Problem 2:")
    t = 60*24*60*60
    t_span = [0, t]
    a = 26553.4 # Semi-major axis [km]
    ecc = 0.741 # Eccentricity
    h = 69084.1 # Angular momentum [km^2/s]
    inc = 63.4 # Inclination [deg]
    raan = 0 # Right ascention of ascending node [deg]
    w = 270 # Arguement of perigee [deg]
    theta = 0 # True anomaly [deg]
    coe = [a, ecc, h, inc, raan, w, theta]
    y0 = coe[1:3] + list(np.deg2rad(coe[3:7]))
    perts = perturbations.Perturbation()
    perts.set_n_body(datetime.datetime.now(), ["Sun"])
    start = time.time()
    sol = solve_ivp(eom.vop, t_span, y0, method='RK45', events=eom.alt_checker_vop, args=(mue, perts,), rtol=1e-8, atol=1e-8)
    elapsed = time.time() - start
    print(f"Time it took to propagate Variation of Parameter: {elapsed} [s]")
    plotter.coe_plot(sol.y[0:6], sol.t/(24*60*60), "coe", mue, re, "Problem 2", deg=False, diff=False)

    print("Over 60 days, there is not a significant change in the orbit. "\
    "The RAAN seems to be the one that is changing the biggest, and it "\
    "seems to be going down without showing any cyclical movement. This "\
    "is concerning because RAAN is the second most expensive correction "\
    "after inclination.")

def problem_3():
    ## Problem 3
    print("Problem 3:")
    mustar = moon.mu/(mue+moon.mu)
    d = {
        "a": {
            "r": np.array([-1.05, 0]),
            "v": np.array([0, -0.066429]),
            "t": 2*np.pi
        },
        "b": {
            "r": np.array([-1.15, 0]),
            "v": np.array([0, -0.0086882909]),
            "t": 29.46
        },
        "c": {
            "r": np.array([0.1, 0]),
            "v": np.array([-3.35, 3]),
            "t": 3.6
        },
        "d": {
            "r": np.array([0.1, 0]),
            "v": np.array([-3.37, 3]),
            "t": 6
        },
        "e": {
            "r": np.array([0.1, 0]),
            "v": np.array([-3.4, 3]),
            "t": 6
        },
        "f": {
            "r": np.array([1.25, 0]),
            "v": np.array([0, 0]),
            "t": 2*np.pi
        },
        "g": {
            "r": np.array([-0.5, 0]),
            "v": np.array([0, 0]),
            "t": 2*np.pi
        },
        "h": {
            "r": np.array([-1.1, 0]),
            "v": np.array([0, 0]),
            "t": 2*np.pi
        }
    }
    # Largrange point
    l_point, mustar = util.lagrange_point(mue, moon.mu)
    l1, l2, l3, l4, l5 = l_point
    # Big primary
    bp = [mustar, 0] 
    # Small primary
    sp = [mustar-1, 0]
    for idx, p in enumerate(d):
        y0 = np.concatenate((d[p]["r"], d[p]["v"]))
        t_span = [0, d[p]["t"]]
        sol = solve_ivp(eom.cr3bp, t_span, y0, method="RK45", args=(mustar,), rtol=1e-8, atol=1e-8)
        for x in range(1,3):
            plt.figure(x)
            plt.subplot(2, 4, idx+1)
            plt.plot(sol.y[0], sol.y[1])
            plt.title(p+": Orbit")
            plt.xlabel("x [DU]")
            plt.ylabel("y [DU]")
            if x == 2:
                plt.plot(bp[0], bp[1], "ko", markersize=10)
                plt.plot(sp[0], sp[1], "ko", markersize=5)
                plt.plot(l1[0], l1[1], "b*", markersize=5)
                plt.plot(l2[0], l2[1], "g*", markersize=5)
                plt.plot(l3[0], l3[1], "c*", markersize=5)
                plt.plot(l4[0], l4[1], "r*", markersize=5)
                plt.plot(l5[0], l5[1], "r*", markersize=5)
                plt.xlim([-1.5, 2])
                plt.ylim([-1.5, 1.5])
    plt.figure(1)
    plt.suptitle("Problem 3 Orbit")
    plt.figure(2)
    plt.suptitle("Problem 3 Lagrange Points and Orbit")
    plt.show()
    print("The trajectories in c.~ e. are called free-return trajectories "\
    "because these trajectories go out of Earth and back into a close "\
    "area near the Earth. Given the correct position and velocity vectors, "\
    "the spacecraft can return to close to Earth orbit without any additional "\
    "delta-V. Plots a. and h. can be used for orbits that want to "\
    "stay near the primary as shown. Plot g. can be used to stay near the bigger primary. "\
    "b. can be used to operate around the big and small primary. f. can be "\
    "used to escape out of the system and away. I decided to plot each "\
    "individually to see what they look like.")

def problem_4():
    ## Problem 4
    print("Problem 4:")
    _ = util.jacobi(earth.mass, moon.mass)
    print("To visualize the L4 and L5 zero velocity curves, I set the "\
    "Jacobi constant to 3. Otherwise, the value would never show up. "\
    "David Eagle also did a similar thing so I assumed this asssumption "\
    "was valid. You can see that the L1 opens up first, then L2, then L3, "\
    "then L4 and L5 at the same time. L4 and L5 are the hardest to get into. "\
    "Lower Jacobi constant means that there is a higher velocity value "\
    "associated to the position so these numbers make sense.")

def main():
    problem_1()
    problem_2()
    problem_3()
    problem_4()

if __name__=="__main__":
    main()
# Homework 2
# AERO 452
# Ryo Takatori

# Modules
import matplotlib.pyplot as plt
import numpy as np
import sympy as sym
from orbit_tools import equation_of_motion as eom
from orbit_tools import planet_data as p_d

earth = p_d.Earth()
mue = earth.mu
re = earth.radius

def problem_1():
    # Problem 1
    print("Problem 1:")
    # 1. At time t=0 a particle is at the origin of a CW frame with a relative
    # velocity of [1 -1 1] (m/s). What will be the relative speed of the
    # particle to the origin after 1/4th of a period?

    # Calculation

    # Initial Relative velocity [m/s]
    delta_v_0 = [1, -1, 1]
    delta_x_dot_0, delta_y_dot_0, delta_z_dot_0 = delta_v_0
    # Since  a = ((((P)/(2*pi))^2)*mu)^(1/3) with P being period, and since we
    # are in CW frame, it will be circular, so a = R. This means that mean
    # motion n, is sqrt(mu/R^3), and therefore, n = 2*pi/P, and t = P/4, so nt
    # is pi/2.
    nt = np.pi/2 # Mean motion*time [rad]

    # Relative velocity [m/s]
    delta_x_dot = delta_x_dot_0*np.cos(nt) + 2*delta_y_dot_0*np.sin(nt)
    delta_y_dot = 4*delta_y_dot_0*np.cos(nt) - 2*delta_x_dot_0*np.sin(nt) - 3*delta_y_dot_0
    delta_z_dot = delta_z_dot_0*np.cos(nt)
    delta_v = np.array([delta_x_dot, delta_y_dot, delta_z_dot])

    # Relative speed [m/s]
    speed = np.linalg.norm(delta_v)

    # Display
    print(f"Relative speed of the particle to the origin after 1/4th of a period: {round(speed, 5)} [m/s]")
    print("I feel good about this answer because it is on the magnitude of a "\
        "few meters per second. It has only been a fourth of a period, so we "\
        "should expect a value like this. I also thought it was interesting "\
        "that if we specify a fraction of a period, we can find the speed "\
        "without knowing much about the orbit, given CW, realtive velocity, and "\
        "inital position of 0.")

def problem_2():
    ## Problem 2
    print("Problem 2:")

    # 2. Curtis 7.15
    # A GEO satellite strikes some orbiting debris and is found 2 h afterward
    # to have drifted to the position dr = -10^i + 10^j [km] relative to its
    # original location. At that time, the only slightly damaged satellite
    # initiates a two-impulse maneuver to return to its original location in 6
    # h. Find the total delta-v for this maneuver.

    # Calculation
    # Period
    T = 24*60*60

    # Mean motion
    n=(2*np.pi)/T

    # Time
    t_drift = 2*60*60 # time to drift
    t_imp = 6*60*60 # time for two impulse
    delta_r_0 = np.array([-10, 10, 0]) # intial relative vector

    # Relative velocity of spacecraft 2h after impact [km/s]
    x = sym.Symbol("x")
    y = sym.Symbol("y")
    eqn1 = (2/n)*(1-np.cos(n*t_drift))*y + (x/n)*np.sin(n*t_drift) + 10
    eqn2 = (2/n)*(np.cos(n*t_drift)-1)*x + y*(-3*t_drift+((4/n)*np.sin(n*t_drift))) - 10
    solution = sym.solve((eqn1, eqn2), (x, y))
    delta_x_dot_0 = np.array(solution[x]).astype(np.float64)
    delta_y_dot_0 = np.array(solution[y]).astype(np.float64)
    delta_x_dot = delta_x_dot_0*np.cos(n*t_drift) + 2*delta_y_dot_0*np.sin(n*t_drift)
    delta_y_dot = 4*delta_y_dot_0*np.cos(n*t_drift) - 2*delta_x_dot_0*np.sin(n*t_drift) - 3*delta_y_dot_0
    delta_z_dot = 0
    delta_v_0 = np.array([delta_x_dot, delta_y_dot, delta_z_dot])
    delta_v, _, _ = eom.two_impluse(delta_r_0, delta_v_0, n, t_imp)
    print(f"Total Delta-V: {delta_v*1000} [m/s]")
    print("Once again, I feel good about this answer because it is in the "\
        "magnitude of a few meters per second. I also feel extra good because "\
        "I checked both the matrix and equation form to verify my answer.")

def problem_3():
    ## Problem 3
    print("Problem 3:")
    # 3. A target is in a 300 km circular geocentric orbit. The chaser is
    # [-100] km behind the target with the same velocity when it initiates a
    # two-impulse maneuver in an effort to rendezvous at the target within 30
    # mins.

    # a) Find the total delta-v.
    print("a.")

    # Calculation
    r = re + 300 # altitude [km]
    n = np.sqrt(mue/r**3) # Mean motion [rad/s]

    # Initial relative position [km]
    r_rel = np.array([-1, 0, 0])

    # Intial relative velocity [km/s]
    v_rel = np.array([0, 0, 0])

    # Two impluse Maneuver for anytime between 30 [min]
    delta_v_list = []
    for t in range(1, 30*60):
        delta_v, _, _ = eom.two_impluse(r_rel, v_rel, n, t)
        delta_v_list.append(delta_v*1000)
    # Finding optimal Delta-V and respective time
    delta_v_best = min(delta_v_list) # Delta-V [m/s]
    t_best = delta_v_list.index(delta_v_best) # Time [min]

    # Display
    print(f"Lowest Delta-V: {delta_v_best} [m/s] with {(t_best+1)/60} [min]")
    plt.plot([t/60 for t in range(1, 30*60)], delta_v_list)
    plt.plot((t_best+1)/60, delta_v_best, "o")
    plt.text(25, 50,"Lowest $\Delta$V")
    plt.title("$\Delta$V Over Time")
    plt.xlabel("Time [min]")
    plt.ylabel("$\Delta$V [m/s]")
    plt.xlim([0,30])
    plt.show()

    print("I feel good about my answer because once again, I see that the "\
        "magnitude of delta-v is in meters per second. I also made sure that "\
        "lowest delta-v was actually happening at the 30 [min] mark because "\
        "the problem statement wanted within 30 [min] which could mean there "\
        "were better times. But the graph shows otherwise, and that the 30 "\
        "[min] is indeed the lowest delta-v.")

def problem_4():
    ## Problem 4
    print("Problem 4:")
    # 4. Spacecraft A and B are in coplanar, circular geocentric orbits. The
    # orbital radii are 8000 km for the target spacecraft, A, and 7000 km for
    # the chaser spacecraft, B. Assume the LVLH frame is with x radially
    # outward and y is in the direction of the velocity of the target. When B
    # is directly below A at perigee, answer these three questions:

    # a) Knowing what you do about the trend analysis, describe the direction
    # of the relative motion.
    print("a.")
    print("Since the target with respect to the chaser is above the chaser, "\
        "the chaser will drift ahead. We also know that the closer and object "\
        "is to the Earth the faster it orbits. This means that the chaser should "\
        "drift ahead. Both statements agree that the chaser will orbit faster "\
        "and drift ahead.")

    # b) Without using the CW equations and the close proximity equation
    # (meaning use the actual relative motion equations), calculate the
    # relative speed of B to A?
    print("b.")

    # Spacecraft A
    r_A = np.array([8000, 0, 0]) # Position [km]
    v_A = np.array([0, np.sqrt(mue/8000), 0]) # Velocity [km/s]
    a_A = np.array([0, -(mue/(8000**2)), 0]) # Acceleration ]km/s^2]

    # Spacecraft B
    r_B = np.array([7000, 0, 0]) # Position [km]
    v_B = np.array([0, np.sqrt(mue/7000), 0]) # Velocity [km/s]
    a_B = np.array([0, 0, -mue/(7000**2)]) # Acceleration ]km/s^2]

    # Five term acceleration calculation
    r_B_A_LVLH ,v_B_A_LVLH, a_B_A_LVLH = eom.fta(r_A, r_B, v_A, v_B, a_A, a_B)
    speed = np.linalg.norm(v_B_A_LVLH) # Speed [km/s]

    print(f"Realtive speed of Spacecraft B to A: {speed} [km/s]")
    print("Knowing that the two circular orbits are coplanar makes this speed "\
        "setup easy because we know exactly where it is without the COE. "\
        "We feel good about this velocity because we expect spacecraft B"\
        "To be orbiting faster because the orbit is closer to Earth. "\
        "Since the problem wanted speed, it is not apparent, but the vector form "\
        "but indeed, spacecraft B has a faster velocity. We also feel good because "\
        "it matched with our trend analysis.")

def main():
    problem_1()
    problem_2()
    problem_3()
    problem_4()

if __name__=="__main__":
    main()
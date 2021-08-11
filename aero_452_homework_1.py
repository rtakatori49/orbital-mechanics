# Homework 1
# AERO 452
# Ryo Takatori

# Modules
import matplotlib.pyplot as plt
import numpy as np
from orbit_tools import classical_orbital_elements as coe
from orbit_tools import coordinate_transformation as c_t
from orbit_tools import direction_cosine_matrix as dcm
from orbit_tools import equation_of_motion as eom
from orbit_tools import planet_data as p_d
from orbit_tools import rendezvous

earth = p_d.Earth()
mue = earth.mu
re = earth.radius

def problem_1():
    # Problem 1
    print("Problem 1:")

    # 1. Given Spacecraft A has these COEs:
    h_A = 51400 # Angular momentum [km^2/s]
    ecc_A = 0.0006387 # Eccentricity
    inc_A = 51.65 # Inclination [deg]
    raan_A = 15 # Right ascension of ascending node [deg]
    w_A = 157 # Argument of perigee [deg]
    theta_A = 15 # True anomaly [deg]

    # And spacecraft B has these COEs:
    h_B = 51398 # Angular momentum [km^2/s]
    ecc_B = 0.0072696 # Eccentricity
    inc_B = 50 # Inclination [deg]
    raan_B = 15 # Right ascension of ascending node [deg]
    w_B = 140 # Argument of perigee [deg]
    theta_B = 15 # True anomaly [deg]

    coe_A = [ecc_A, h_A, inc_A, raan_A, w_A, theta_A]
    r_A, v_A = coe.coe_to_state(mue, coe_A)
    a_A = h_A**2/((1-ecc_A**2)*mue) # Semi-major axis [km]
    T_A = 2*np.pi*np.sqrt(a_A**3/mue) # Period [s]

    coe_B = [ecc_B, h_B, inc_B, raan_B, w_B, theta_B]
    r_B, v_B = coe.coe_to_state(mue, coe_B)
    a_B = h_B**2/((1-ecc_B**2)*mue) # Semi-major axis [km]
    T_B = 2*np.pi*np.sqrt(a_B**3/mue) # Period [s]

    # a) Calculate the relative position of S/C B to S/C A for 10 periods of S/C A
    print("a.")
    r_B_A_lvlh_list = [[], [], []]
    r_B_A_list = [[], [], []]
    r_A_list = [[], [], []]
    r_B_list = [[], [], []]
    norm_r_B_A_list = []
    time_list = []
    for x in range(int(T_A)):
        r_A, v_A = eom.kepler(r_A, v_A, 10, mue)
        r_B, v_B = eom.kepler(r_B, v_B, 10, mue)
        r_B_A = r_B - r_A
        norm_r_B_A_list.append(np.linalg.norm(r_B_A))
        _, q = c_t.eci_to_lvlh(r_A, v_A)
        r_B_A_lvlh = np.dot(q, r_B_A)
        time_list.append(x*10/(60*60))
        for y in range(3):
            r_A_list[y].append(r_A[y])
            r_B_list[y].append(r_B[y])
            r_B_A_list[y].append(r_B_A[y])
            r_B_A_lvlh_list[y].append(r_B_A_lvlh[y])

    # b) What is the closest approach of B to A and at what time does this occur?
    print("b.")

    r_B_A_close = min(norm_r_B_A_list) # Closest distance [km]
    t_B_A_close = norm_r_B_A_list.index(r_B_A_close)*10 # Closest time [s]
    print(f"The closest approach of B to A: {round(r_B_A_close,3)} [km]")
    print(f"The closest approach time of B to A: {round(t_B_A_close/(60*60),3)} [hr]")
    print("Since we are not applying any kind of change to the spacecraft in "\
        "orbit, we do not expect to see a huge conversion in the distance of "\
        "approach. It is also interesting that the closest approach happens "\
        "very close to the end of 10 periods.")

    # c) Graph the results
    print("c.")

    fig = plt.figure()

    ax = fig.add_subplot(2, 2, 1, projection="3d")
    plt.plot(r_A_list[0], r_A_list[1], r_A_list[2])
    plt.plot(r_B_list[0], r_B_list[1], r_B_list[2])
    plt.title("Orbit of Spacecraft in ECI for 10 Periods")
    plt.xlabel("x [km]")
    plt.ylabel("y [km]")
    ax.set_zlabel("z [km]")
    plt.legend(["Spacecraft A", "Spacecraft B"])

    ax = fig.add_subplot(2, 2, 2, projection="3d")
    plt.plot(r_B_A_list[0], r_B_A_list[1], r_B_A_list[2])
    plt.title("Relative Position of Spacecraft B relative to Spacecraft A in ECI for 10 Periods")
    plt.xlabel("x [km]")
    plt.ylabel("y [km]")
    ax.set_zlabel("z [km]")

    ax = fig.add_subplot(2, 2, 3, projection="3d")
    plt.plot(r_B_A_lvlh_list[0], r_B_A_lvlh_list[1], r_B_A_lvlh_list[2])
    plt.title("Relative Position of Spacecraft B relative to Spacecraft A in LVLH for 10 Periods")
    plt.xlabel("x [km]")
    plt.ylabel("y [km]")
    ax.set_zlabel("z [km]")

    plt.subplot(2,2,4)
    plt.plot(time_list, norm_r_B_A_list)
    plt.title("Distance of Spacecraft B relative to Spacecraft A for 10 periods")
    plt.xlabel("Time [hr]")
    plt.ylabel("Distance [km]")

    plt.suptitle("Problem 1")
    plt.show()

    print("As noticed by the graph, the orbits for the two spacecraft are "\
        "very similar. This is also reflected in the closest approach graph. "\
        "The two spacecraft do get gradually closer over time, but it is not "\
        "very close. It is also continually getting closer, so perhaps if we "\
        "run the calculation for a longer time, we will see it get even closer.")

def problem_2():
    ## Problem 2
    print("Problem 2:")
    # 2. Curtis problem 7.1
    # Two manned spacecraft, A and B (see the figure), are in circular polar (i
    # = 90 [deg]) orbits around the earth. A"s orbital altitude is 300 [km]
    # B"s is 250 [km]. At the instant shown (A over the equator, B over the
    # North Pole), calculate (a) the position, (b) velocity, and (c) the
    # acceleration of B relative to A. A"s y-axis points always in the flight
    # direction, and its x-axis is directed radially outward at all times.

    # Spacecraft A
    inc_A = 90 # Inclination [deg]
    ecc_A = 0 # Eccentricity
    a_A = 300 + re # Semi-major axis [km]
    r_A = np.array([0, a_A, 0]) # Position [km]
    v_A = np.array([0, 0, np.sqrt(mue/a_A)]) # Velocity [km/s]
    acc_A = np.array([0, -(mue/(a_A**2)), 0]) # Acceleration ]km/s^2]

    # Spacecraft B
    inc_B = 90 # Inclination [deg]
    ecc_B = 0 # Eccentricity
    a_B = 250 + re # Semi-major axis [km]
    r_B = np.array([0, 0, a_B]) # Position [km]
    v_B = np.array([0, -np.sqrt(mue/a_B), 0]) # Velocity [km/s]
    acc_B = np.array([0, 0,-mue/(a_B**2)])  # Acceleration ]km/s^2]

    # a.
    print("a.")
    # Calculation
    # Relative position
    r_B_A = r_B - r_A # Position [km]
    norm_r_A = np.linalg.norm(r_A) # Position magnitude [km]
    h_A = np.cross(r_A,v_A) # Angular momentum [km^2/s]
    _, q = c_t.eci_to_lvlh(r_A, v_A)
    r_B_A_LVLH = np.dot(q, r_B_A) # Position in LVLH [km]

    #Display
    print(f"Relative position of Spacecraft B relative to Spacecraft A: {np.around(r_B_A_LVLH, decimals=3)} [km]")
    print("We feel good about this answer, since by inspection of the drawing "\
        "we know that this vector is correct.")

    # b.
    print("b.")
    # Calculation
    # Relative velocity
    omega = h_A/norm_r_A**2 # Relative angular velocity [rad/s]
    v_B_A = v_B - v_A - np.cross(omega,r_B_A) # Relative velocity [km/s]
    v_B_A_LVLH = np.dot(q, v_B_A) # Relative velocity in LVLH [km/s]

    #Display
    print(f"Relative velocity of Spacecraft B relative to Spacecraft A: {v_B_A_LVLH} [km/s]")
    print("We feel good about this answer, since the two objects are in the "\
        "same orbit with relatively close radius, so we expect to see this "\
        "small value in the relative velocity.")

    # c.
    print("c.")
    # Calculation
    # Relative acceleration
    omega_dot_2 = ((-2*np.dot(v_A,r_A))*omega)/norm_r_A**2 # Relative angular acceleration [rad/s^2]
    a_B_A = acc_B - acc_A - np.cross(omega_dot_2, r_B_A) - np.cross(omega,np.cross(omega,r_B_A)) - 2*np.cross(omega,v_B_A) # Relative acceleration [km/s^2]
    a_B_A_LVLH = np.dot(q, a_B_A) # Relative acceleration in LVLH [km/s^2]

    # Display
    print(f"Relative acceleration of Spacecraft B relative to Spacecraft A: {a_B_A_LVLH} [km/s^2]")
    print("We feel good about this answer because the acceleration of an "\
        "orbit is dependant on the semi-major axis, and in this case, the "\
        "circular orbits do not have a big difference in the radius, so we "\
        "expect to see this kind of minute difference in the acceleration.")

def problem_3():
    ## Problem 3
    print("Problem 3:")
    # 3. Chief satellite is at:
    rp = 250 + re # Perigee [km]
    ecc_A = 0.1 # Eccentricity
    a_A = rp/(1-ecc_A) # Semi-major axis [km]
    h_A = np.sqrt(a_A*((1-ecc_A**2)*mue)) # Angular momentum [km^2/s]
    inc_A = 51 # Inclination [deg]
    raan_A = 0 # Right ascension of ascending node [deg]
    w_A = 0 # argument of perigee [deg]
    theta_A = 0 # true anomaly [deg]
    T_A = 2*np.pi*np.sqrt(a_A**3/mue) # Period [s]
    coe_A = [ecc_A, h_A, inc_A, raan_A, w_A, theta_A]
    r_A, v_A = coe.coe_to_state(mue, coe_A)
    # Deputy is at:
    r_rel = np.array([-1, -1, 0]) # Relative position [km]
    delta_v_1 = np.array([0, 2, 0])/1000 #Relative velocity [km/s]

    t_span = [0, T_A*10] # timespan [s]
    y0 = np.concatenate((r_A, v_A, r_rel, delta_v_1)) # initial states
    rendezvous.eom_leom(t_span, y0, np.array([0, 0, 0]), False, False, 0, "Relative separation of two object over 10 periods", True)
    

    print("Seeing the graphs of the relative position for the x and y "\
    "direction, it seems that the separation is getting larger and larger "\
    "as time progresses. However, the separation does not blow up, "\
    "because the relative velocity difference was very small, so we expect "\
    "to see this kind of slow separation.")

def problem_4():
    ## Problem 4
    print("Problem 4:")
    # 4. Curtis 7.7
    # A space station is in a 90-min period earth orbit. At t = 0, a satellite
    # has the following position and velocity components relative to a CW frame
    # attached to the space station: dr =^i [km], dv = 10^j [m/s]. How far is
    # the satellite from the space station 15 [min] later?

    # Calculation
    # Position of spacecraft with respect to space station [m]
    delta_x_0 = 1
    delta_y_0 = 0
    delta_z_0 = 0
    delta_r_0 = np.array([delta_x_0, delta_y_0, delta_z_0])

    # Velocity of spacecraft with respect to space station [m/s]
    delta_x_dot_0 = 0
    delta_y_dot_0 = 10/1000
    delta_z_dot_0 = 0
    delta_v_0 = np.array([delta_x_dot_0, delta_y_dot_0, delta_z_dot_0])

    # Period
    T = 90*60

    # Mean motion
    n=(2*np.pi)/T

    # Time
    t = 15*60

    # CW matrix
    phi_rr, phi_rv, phi_vr, phi_vv = dcm.clohessy_wiltshire(n, t)

    # Distance after time [km]
    delta_r = np.dot(phi_rr, delta_r_0) + np.dot(phi_rv, delta_v_0)
    dist = np.linalg.norm(delta_r)
    delta_v = np.dot(phi_vr, delta_r_0) + np.dot(phi_vv, delta_v_0)
    speed = np.linalg.norm(delta_v)

    # Display
    print(f"Statellite from the space station 15 [min] later: {round(dist, 5)} [km]")
    print(f"Statellite velocity from the space station 15 [min] later: {round(speed, 5)} [km/s]")
    print("Since we know that the spacecraft started 1000 [km] away, "\
        "and the resulting separation distance is 11.22 [km], it is safe "\
        "to say that the spacecraft is approaching the spacestation. However "\
        "I assume that this is not an intended approach, as closing a gap of "\
        "roughly 990 [km] in 15 [min] seems way too fast for a dock.")

def main():
    problem_1()
    problem_2()
    problem_3()
    problem_4()

if __name__=="__main__":
    main()
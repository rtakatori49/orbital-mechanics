# Homework 1
# AERO 557
# Ryo Takatori

# Modules
import datetime
from matplotlib.pyplot import plot
import numpy as np
from orbit_tools import classical_orbital_elements as coe
from orbit_tools import orbit_determination as od
from orbit_tools import lambert
from orbit_tools import planet_data as p_d
from orbit_tools import plotter

earth = p_d.Earth()
mue = earth.mu

def problem_1():
    # Problem 1
    print("Problem 1:")
    # Given
    # Observation station
    lat = 40 # Latitude [deg]
    long = -110 # Longitude [deg]
    alt = 2 # Altitiude [km]

    # Observation 1
    t_1 = datetime.datetime(2010, 8, 20, 11, 30, 0)
    ra_1 = -33.0588410 # Right Ascension [deg]
    dec_1 = -7.2056382 # Declination [deg]

    # Observation 2
    t_2 = datetime.datetime(2010, 8, 20, 11, 50, 0)
    ra_2 = 55.0931551 # Right Ascension [deg]
    dec_2 = 36.5731946 # Declination [deg]

    # Observation 3
    t_3 = datetime.datetime(2010, 8, 20, 12, 0, 0)
    ra_3 = 98.7739537 # Right Ascension [deg]
    dec_3 = 31.1314513 # Declination [deg]

    # All observation combined for easier data handling
    datetime_list = [t_1, t_2, t_3]
    ra_list = [ra_1, ra_2, ra_3]
    dec_list = [dec_1, dec_2, dec_3]
    
    # Gauss" Method (non-extended)
    r2, v2 = od.gauss(lat, long, alt, datetime_list, ra_list, dec_list, extend=False)
    norm_r2 = np.linalg.norm(r2)
    norm_v2 = np.linalg.norm(v2)

    # Display
    print("Gauss'"" Method (non-extended):")
    print(f"r2 Vector: {r2} [km]")
    print(f"r2 Magnitude: {norm_r2} [km]")
    print(f"v2 Vector: {v2} [km/s]")
    print(f"v2 Magnitude: {norm_v2} [km/s]")
    print("Classical Orbital Elements:")
    coe_data = coe.state_to_coe(mue, r2, v2)
    plotter.coe_table(coe_data)
    
    # Gauss" Method (extended)
    r2_e, v2_e = od.gauss(lat, long, alt, datetime_list, ra_list, dec_list, extend=True)
    norm_r2_e = np.linalg.norm(r2_e)
    norm_v2_e = np.linalg.norm(v2_e)
    
    # Display
    print("Gauss'"" Method (extended):")
    print(f"r2 Vector: {r2_e} [km]")
    print(f"r2 Magnitude: {norm_r2_e} [km]")
    print(f"v2 Vector: {v2_e} [km/s]")
    print(f"v2 Magnitude: {norm_v2_e} [km/s]")
    print("Classical Orbital Elements:")
    coe_data_e = coe.state_to_coe(mue, r2_e, v2_e)
    plotter.coe_table(coe_data_e)

    # Double-R Method
    r2_dr, v2_dr = od.double_r(lat, long, alt, datetime_list, ra_list, dec_list)
    norm_r2_dr = np.linalg.norm(r2_dr)
    norm_v2_dr = np.linalg.norm(v2_dr)
    
    # Display
    print("Double-R Method:")
    print(f"r2 Vector: {r2_dr} [km]")
    print(f"r2 Magnitude: {norm_r2_dr} [km]")
    print(f"v2 Vector: {v2_dr} [km/s]")
    print(f"v2 Magnitude: {norm_v2_dr} [km/s]")
    print("Classical Orbital Elements:")
    coe_data_dr = coe.state_to_coe(mue, r2_dr, v2_dr)
    plotter.coe_table(coe_data_dr)

    print("We can see that the Gauss"" extended version and double-r " \
    "iteration seems to converge to the same value. This shows that " \
    "these two methods are much better at inital orbit determination " \
    "than the non-extended Gauss. This is expected as the Gauss"" non-" \
    "extended does not have any form of iteration which means it does " \
    "not refine the answer. The Gauss"" non-extended was off by more " \
    "than 300 [km] which is quite significant, and furthermore, the " \
    "velocity value was off by 0.4 [km/s] which is also siginificantly " \
    "off. These differences show up huge in the orbital elements. " \
    "the values between Gauss extended and double-r are not noticable. " \
    "However, in the non-extended, it is clear that the orbit found was " \
    "different from the one found in the other two methods.")

def problem_2():
    # Problem 2
    print("Problem 2:")

    # Given
    # Position [km]
    r0 = np.array([15945.34, 0, 0]) # 1
    rf = np.array([12214.83899, 10249.46731, 0]) # 2
    # Short way
    tm = 1
    # Time of transfer [min]
    dt = 76.0
    # Number of revs: 1
    m = 0

    # Minimum Energy
    _, _, dtmin_me, _, v0_me = lambert.minimum_energy(r0, rf, tm, mue)
    print("Minimum Energy:")
    print(f"v0 Vector: {v0_me} [km/s]")
    print(f"Minimum energy transfer time: {dtmin_me/60} [min]")

    # Gauss
    v0_lg, vf_lg = lambert.gauss(r0, rf, dt*60, tm, mue)
    print("Gauss:")
    print(f"v0 Vector: {v0_lg} [km/s]")
    print(f"v1 Vector: {vf_lg} [km/s]")

    # Universal Variables
    v0_luv, vf_luv = lambert.universal_variable(r0, rf, dt*60, mue)
    print("Universal Method:")
    print(f"v0 Vector: {v0_luv} [km/s]")
    print(f"v1 Vector: {vf_luv} [km/s]")

    # # Izzo/Gooding Method
    # v0_ig, vf_ig, extremal_distances_ig_2, exitflag_ig_2] = ...
    # lambert.iz(r0_2, r1_2, dt_2/(60*24), m_2, mue);
    # print("Izzo/Gooding Method:")
    # print(f"v0 Vector: {v0_ig} [km/s]")
    # print(f"v1 Vector: {vf_ig} [km/s]")

    print("All four methods show similar results, and particularly the " \
        "universal method and Izzo/Gooding method. We could not take " \
        "advantage of the Izzo/Gooding method in this problem since we were " \
        "asked only about a single revolution. Izzo/Gooding method would have " \
        "been far superior in the case of 2+ revolutions. In this specific" \
        "problem, it seems like the time of transfer found is close to the " \
        "value given in the problem, but this is by coincidence and if we " \
        "were to conduct this with some other value this would not be the case.")


def main():
    problem_1()
    problem_2()

if __name__=="__main__":
    main()
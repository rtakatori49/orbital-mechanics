import datetime
import pandas as pd
from orbit_tools import equation_of_motion as eom
from orbit_tools import orbit_determination as od
from orbit_tools import planet_data as p_d

earth = p_d.Earth()
mue = earth.mu
def read_text(filename):
    with open(filename) as f:
        raw_data = f.readlines()
    datetime_list = []
    range_list = []
    azimuth_list = []
    elevation_list = []
    for line in raw_data:
        data = line.split()
        second_string = data[5].split(".")
        datetime_list.append(datetime.datetime(int(data[0]), int(data[1]),
            int(data[2]), int(data[3]), int(data[4]), int(second_string[0]),
            int(second_string[1])))
        azimuth_list.append(float(data[6]))
        elevation_list.append(float(data[7]))
        range_list.append(float(data[8]))
    return datetime_list, range_list, azimuth_list, elevation_list

def problem_1():
    # Problem 1
    print("Problem 1:")
    # Original samples
    print("Original samples:")
    x_o_i_original = range(1, 11)
    y_o_i_original = [1, 2, 2, 3, 4, 4, 5, 7, 9, 9]
    od.lls(x_o_i_original, y_o_i_original)
    print("Altered samples:")
    x_o_i_altered = [1, 2, 3, 4, 5, 6, 7, 8, 10]
    y_o_i_altered = [1, 2, 2, 3, 4, 4, 5, 7, 9]
    od.lls(x_o_i_altered, y_o_i_altered)
    print("Although the highest resdiual was lower than two times the "\
    "sigma, I pulled the highest residual value from the sample and "\
    "ran the LLS with the fewer sample size. As was discussed in class "\
    "the root mean square saw a lower value, but the confidence interval "\
    "also increased. This means that we have a overall better fitting "\
    "set of data, but we have a lower confidence in those values due to "\
    "decrease in the sample size we have. Often times it is better to just "\
    "increase the sample size rather than getting rid of one bad value.")

def problem_2():
    # Problem 2
    print("Problem 2:")
    # datetime_list, range_list, az_list, el_list\
    #     = read_text("Data4Problem2Hw2_2020test.txt")
    df_obs = pd.read_csv("Data4Problem2Hw2_2020test.csv")
    lat = 21.5748 # Latitude [deg]
    long = -158.2706 # Longitude [deg]
    alt = 300.2/1000 # Altitude [m]
    rho_error = 92.5/1000 # Range error [m]
    az_error = 0.0224 # Azimuth error [deg]
    el_error = 0.0139 # Elevation error [deg]
    od.nlwls(lat, long, alt, df_obs, rho_error, az_error, el_error, eom.two_body, mue)


def problem_3():
    # Problem 3
    print("Problem 3:")

def main():
    #problem_1()
    problem_2()
    problem_3()

if __name__=="__main__":
    main()
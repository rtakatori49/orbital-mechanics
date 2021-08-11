import datetime
from orbit_tools import planet_data as p_d
from orbit_tools import utility as util

def problem_1():
    # Problem 1
    print("Problem 1:")
    dep_early = datetime.datetime(2005, 6, 6, 0, 0, 0) # Earth departure earliest date (June 6th 2005)
    dep_late = datetime.datetime(2005, 11, 7, 0, 0, 0) # Earth departure latest date (November 7th 2005)
    arr_early = datetime.datetime(2005, 12, 1, 0, 0, 0) # Mars arrival earliest date (December 1st 2005)
    arr_late = datetime.datetime(2007, 2, 24, 0, 0, 0) # Mars arrival latest date (February 24th 2007)
    _, _, _ = util.porkchop(p_d.Earth(), dep_early, dep_late,
        p_d.Mars(), arr_early, arr_late, 5)

def problem_2():
    # Problem 2
    print("Problem 2:")
    print("Word problem done on paper.")

def problem_3():
    # Problem 3
    print("Problem 3:")

def main():
    problem_1()
    problem_2()
    problem_3()

if __name__=="__main__":
    main()
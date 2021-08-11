# Utility
# Functions useful or necessary for orbital mechanics

# Modules
import datetime
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from . import classical_orbital_elements as coe
from . import equation_of_motion as eom
from . import lambert
from . import planet_data as p_d
from . import time_trans as t_t
from . import misc_data
au = misc_data.au
earth = p_d.Earth()
mue = earth.mu # Earth gravitational constant [km^3/s^2]
re = earth.radius # Earth radius [km]


# Two-Line Element Reader
# Ryo Takatori
# 11/17/2020
def tle_reader(filename):
    # Open TLE text file
    with open(filename, 'r') as f:
        tle_lines_raw = f.readlines()
    # Preallocate list
    tle_title = []
    tle_1 = []
    tle_2 = []
    tle = []

    # Sort TLE
    tle_title.append(" ".join(tle_lines_raw[0].split())) # title
    tle.append(tle_title)
    # Line 1
    tle_1.append(int(tle_lines_raw[1][0])) # line number
    tle_1.append(int(tle_lines_raw[1][2:7])) # satellite catalog number
    tle_1.append(tle_lines_raw[1][7]) # classification
    # International designator
    tle_1.append(int(tle_lines_raw[1][9:11])) # last two digits of launch year
    tle_1.append(int(tle_lines_raw[1][11:14])) # launch number of the year
    tle_1.append(" ".join(tle_lines_raw[1][14:17].split())) # piece of launch
    #
    tle_1.append(int(tle_lines_raw[1][18:20])) # epoch year (last two digits of launch year)
    tle_1.append(float(tle_lines_raw[1][20:32])) # epoch (day of the year and fractional portion of the day)
    tle_1.append(float(tle_lines_raw[1][33:43])) # first derivative of mean motion (ballistic coefficient)
    tle_1.append(float(tle_lines_raw[1][43]+"0."+" ".join(tle_lines_raw[1][44:50].split()))*10**int(tle_lines_raw[1][50:52])) # second derivative of mean motion
    tle_1.append(float(tle_lines_raw[1][53]+"0."+" ".join(tle_lines_raw[1][54:59].split()))*10**int(tle_lines_raw[1][59:61])) # BSTAR
    tle_1.append(int(tle_lines_raw[1][62])) # ephemeris type
    tle_1.append(int(tle_lines_raw[1][64:68])) # element set number
    tle_1.append(int(tle_lines_raw[1][68])) # check sum
    tle.append(tle_1)
    # Line 2
    tle_2.append(int(tle_lines_raw[2][0])) # line number
    tle_2.append(int(tle_lines_raw[2][2:7])) # satellite catalog number
    tle_2.append(float(tle_lines_raw[2][8:16])) # inclination [deg]
    tle_2.append(float(tle_lines_raw[2][17:25])) # right ascension of ascending node [deg]
    tle_2.append(float("0."+tle_lines_raw[2][26:33])) # eccentricity
    tle_2.append(float(tle_lines_raw[2][34:42])) # argument of perigee [deg]
    tle_2.append(float(tle_lines_raw[2][43:51])) # mean anomaly [deg]
    tle_2.append(float(tle_lines_raw[2][52:63])) # mean motion [rev/day]
    tle_2.append(int(tle_lines_raw[2][63:68])) # revolution number at epoch [rev]
    tle_2.append(int(tle_lines_raw[2][68])) # check sum
    tle.append(tle_2)
    return tle

# Class for Satellite
# Ryo Takatori
# 11/17/2020
class Satellite:
    def __init__(self, filename):
        # Two-Line Element Set
        self.tle = tle_reader(filename)
        self.name = self.tle[0][0] # name
        # Line 1
        self.line_number_1 = self.tle[1][0] # line number
        self.catalog_number_1 = self.tle[1][1] # satellite catalog number
        self.classification = self.tle[1][2] # classification
        # International designator
        self.id_launch_year = self.tle[1][3] # last two digits of launch year
        self.id_launch_number = self.tle[1][4] # launch number of the year
        self.id_launch_piece = self.tle[1][5] # piece of launch
        #
        self.epoch_year = self.tle[1][6] # epoch year (last two digits of launch year)
        self.epoch_day_of_year = self.tle[1][7] # epoch (day of the year and fractional portion of the day)
        self.ballistic_coefficient = self.tle[1][8] # first derivative of mean motion (ballistic coefficient)
        self.d2_mean_motion = self.tle[1][9] # second derivative of mean motion
        self.bstar = self.tle[1][10] # BSTAR
        self.ephemeris_type = self.tle[1][11] # ephemeris type
        self.element_set_number = self.tle[1][12] # element set number
        self.checksum_1 = self.tle[1][13] # check sum
        # Line 2
        self.line_number_2 = self.tle[2][0] # line number
        self.catalog_number_2 = self.tle[2][1] # satellite catalog number
        self.inclination = self.tle[2][2] # inclination [deg]
        self.right_ascension_ascending_node = self.tle[2][3] # right ascension of ascending node [deg]
        self.eccentricity = self.tle[2][4] # eccentricity
        self.argument_of_perigee = self.tle[2][5] # argument of perigee [deg]
        self.mean_anomaly = self.tle[2][6] # mean anomaly [deg]
        self.mean_motion = self.tle[2][7] # mean motion [rev/day]
        self.revolution = self.tle[2][8] # revolution number at epoch [rev]
        self.checksum_2 = self.tle[2][9] # check sum

        # Classical Orbital Elements and other related elements and parameters
        self.n = (self.mean_motion*2*np.pi)/(24*60*60) # mean motion
        self.me = self.mean_anomaly # mean anomaly [deg]
        self.a = (mue/self.n**2)**(1/3) # semi-major axis [km]
        self.T = 2*np.pi*np.sqrt((self.a**3)/mue) # period [s]
        self.ecc = self.eccentricity # eccentricity
        self.rp = self.a*(1-self.ecc) # perigee [km]
        self.ra = (2*self.a) - self.rp # apogee [km]
        self.zp = self.rp - re # perigee altitude [km]
        self.za = self.ra - re # apogee altitude [km]
        self.E = eom.kepler_e(self.ecc, np.deg2rad(self.me)) # eccentric anomaly
        self.h = np.sqrt(self.a*mue*(1-self.ecc**2)) # specifc engular momentum [km^2/s]
        self.inc = self.inclination # inclination [deg]
        self.raan = self.right_ascension_ascending_node # right ascension of ascending node [deg]
        self.w = self.argument_of_perigee # argument of perigee [deg]
        self.theta = np.rad2deg(2*np.arctan(np.sqrt(((1+self.ecc)/(1-self.ecc)))*np.tan(np.rad2deg(self.E)/2))) # true anomaly [deg]
        self.coe = [self.a, self.ecc, self.h, self.inc, self.raan, self.w, self.theta]
        # Vector
        r, v = coe.coe_to_state(mue, self.coe[1:7]) # coe to vectors
        self.r = r # position vector [km]
        self.v = v # velocity vector [km/s]
        # Check for past 2000 or before
        if self.epoch_year > 56:
            epoch_full_year = 1900 + self.epoch_year
        else:
            epoch_full_year = 2000 + self.epoch_year
        self.epoch = datetime.datetime(epoch_full_year, 1, 1) + datetime.timedelta(self.epoch_day_of_year - 1) # epoch in standard datetime

# Planetary Ephemerides from Meeus (1991:202-204) and J2000.0
def planetary_elements(planet, t):
    planet_coe = planet.ephemeris(t)
    a, ecc, inc, raan, w_hat, L = planet_coe
    # Convert to km:
    sun = p_d.Sun()
    # Calculate specific h, ta, w
    h = np.sqrt(a*au*sun.mu*(1-ecc**2)) # Specific angular momentum [km^2/s]
    M = np.deg2rad(L-w_hat) # Mean anomaly [radian]
    E = eom.kepler_e(ecc, M) # Eccentric anomaly [radian]
    theta = 2*np.arctan(np.sqrt((1+ecc)/(1-ecc))*np.tan(E/2)) # True anomaly [degree]
    theta = np.rad2deg(theta)
    if theta < 0:
        theta = 360 + theta
    w = w_hat - raan # Argument of perigee [degree]
    coe_data = [ecc, h, inc, raan, w, theta]
    r, v = coe.coe_to_state(sun.mu, coe_data)
    return planet_coe, coe_data, r, v

# Porkchop Plot
def porkchop(dep_planet, dep_early, dep_late,
    arr_planet, arr_early, arr_late, step):
    sun = p_d.Sun()
    d_e = t_t.jd(dep_early)
    d_l = t_t.jd(dep_late)
    a_e = t_t.jd(arr_early)
    a_l = t_t.jd(arr_late)
    dep_list = np.arange(d_e, d_l, step)
    arr_list = np.arange(a_e, a_l, step)
    dt = np.zeros((len(arr_list), len(dep_list)))
    vinf = np.zeros((len(arr_list), len(dep_list)))
    c3 = np.zeros((len(arr_list), len(dep_list)))
    for x, arr in enumerate(arr_list):
        _, _, r_arr, v_arr = planetary_elements(arr_planet, arr)
        for y, dep in enumerate(dep_list):
            _, _, r_dep, v_dep = planetary_elements(dep_planet, dep)
            dt_temp = arr - dep
            dt[x, y] = dt_temp
            v_dep_trans, v_arr_trans = lambert.universal_variable(r_dep, r_arr,
                dt_temp*24*60*60, sun.mu)
            vinf[x, y] = np.linalg.norm(v_arr-v_arr_trans)
            c3[x, y] = np.linalg.norm(v_dep_trans-v_dep)**2
    # Display
    x = np.arange(0, d_l-d_e, step)
    y = np.arange(0, a_l-a_e, step)
    X, Y = np.meshgrid(x, y)
    levelc3 = [round(level, 1) for level in np.logspace(0, 2, 20)]
    levelvinf = [round(level, 1) for level in np.logspace(0, 1, 10)]
    leveldt = np.linspace(50, 500, 10)
    C3 = plt.contour(X, Y, c3, levels=levelc3, colors='r')
    Vinf = plt.contour(X, Y, vinf, levels=levelvinf, colors='b')
    DT = plt.contour(X, Y, dt, levels=leveldt, colors='k')
    plt.clabel(C3, inline=True, fontsize=10)
    plt.clabel(Vinf, inline=True, fontsize=10)
    plt.clabel(DT, inline=True, fontsize=10)
    plt.title(f"Pork Chop Plot ({dep_planet.name} => {arr_planet.name}) ")
    plt.xlabel(f"Days Past {dep_early}")
    plt.ylabel(f"Days Past {arr_early}")
    plt.title(f"Pork Chop Plot ({dep_planet.name} => {arr_planet.name}) ")
    plt.xlabel(f"Days Past {dep_early}")
    plt.ylabel(f"Days Past {arr_early}")
    plt.legend([C3.collections[0], Vinf.collections[0], DT.collections[0]],
        ["$C_3 [km^2/s^2]$", "$V_\infty @ Mars [km/s]$",
        "Time of Flight [days]"], loc="upper left")
    plt.show()
    return c3, vinf, dt

# Lagrange Point (Left Hand Side)
def lagrange_point(mb, ms):
    # mu^* calculation
    # mb: mass of bigger primary
    # ms: mass of smaller primary
    mustar = ms/(mb+ms)

    def point(eqn, symbol):
        solution = sym.solve(eqn, symbol)
        l = [sol.real for sol in [complex(sol) for sol in solution]
            if sol.imag == 0]
        return [l[0], 0]

    # L1
    l1point = sym.Symbol("l1point")
    eqn1 = l1point \
        + ((1-mustar)/((l1point-mustar)**2)) \
        - (mustar/((l1point-mustar+1)**2))
    l1 = point(eqn1, l1point)

    # L2
    l2point = sym.Symbol("l2point")
    eqn2 = l2point \
        + ((1-mustar)/((l2point-mustar)**2)) \
        + (mustar/((l2point-mustar+1)**2))
    l2 = point(eqn2, l2point)

    # L3
    l3point = sym.Symbol("l3point")
    eqn3 = l3point \
        - ((1-mustar)/((l3point-mustar)**2)) \
        - (mustar/((l3point-mustar+1)**2))
    l3 = point(eqn3, l3point)

    # L4 & L5
    x = mustar - 0.5
    y = np.sqrt(3)/2
    l4 = [x, y]
    l5 = [x, -y]
    return [l1, l2, l3, l4, l5], mustar

## Jacobi Constant (Left Hand Side)
def jacobi(mb, ms):
    # Lagrange point
    l_point, mustar = lagrange_point(mb,ms)
    l1, l2, l3, l4, l5 = l_point
    # Pre-allocation
    limit = []
    # Jacobi constant at Lagrange points
    for idx, point in enumerate(l_point):
        lx = point[0]
        ly = point[1]
        r1 = np.sqrt(((lx-mustar)**2)+(ly**2))
        r2 = np.sqrt(((lx+1-mustar)**2)+(ly**2))
        if idx > 2:
            temp = 3
        else:
            temp = lx**2 + ly**2 + ((2*(1-mustar))/r1) + ((2*mustar)/r2)
        limit.append(temp)

    # Big primary
    bp = [mustar, 0]
    # Small primary
    sp = [mustar-1 ,0]
    # Jacobi Constant
    x = y = np.arange(-2, 2, 0.01)
    X, Y = np.meshgrid(x,y)
    r1 = np.sqrt(((X-mustar)**2)+(Y**2))
    r2 = np.sqrt(((X+1-mustar)**2)+(Y**2))
    C = X**2+Y**2+((2*(1-mustar))/r1) + ((2*mustar)/r2)

    # Display
    plt.plot(bp[0], bp[1], "ko", markersize=10)
    plt.plot(sp[0], sp[1], "ko", markersize=5)
    plt.plot(l1[0], l1[1], "b*", markersize=5)
    plt.plot(l2[0], l2[1], "g*", markersize=5)
    plt.plot(l3[0], l3[1], "c*", markersize=5)
    plt.plot(l4[0], l4[1], "r*", markersize=5)
    plt.plot(l5[0], l5[1], "r*", markersize=5)
    plt.xlim([-2, 2])
    plt.ylim([-2, 2])
    CS = plt.contour(X, Y, C, levels=np.linspace(3, 4, 5))
    plt.clabel(CS, inline=True, fontsize=10)
    plt.title("C values")
    plt.xlabel("x [DU]")
    plt.ylabel("y [DU]")
    plt.legend(["Big Primary", "Small Primary", 
        "$L_1$", "$L_2$", "$L_3$", "$L_4$", "$L_5$"])
    plt.show()
    return C
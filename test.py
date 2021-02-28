# Test script to test all functions

import orbit_tools.dcm as dcm
import orbit_tools.coe as coe
import orbit_tools.time_trans as time_trans
import orbit_tools.celestial_body_ephemeris as cbe
import perturbations
import numpy as np
import datetime
import orbit_tools.tle_reader as tle_reader
import orbit_tools.satellite as satellite
import orbit_tools.planet_data as planet_data
from scipy.integrate import solve_ivp
import orbit_tools.eom as eom
import orbit_tools.plotter as plotter
cx = dcm.rotx(np.pi/4)
print(cx)

r = [-10766.25, -3383.89, 9095.35]
v = [-0.9250, -5.1864, -2.1358]
earth_mu = 398600
earth_radius = 6378
# r_mag = earth_radius + 100.0
# v_mag = np.sqrt(earth_mu/r_mag)
# r = [r_mag, 0, 0]
# v = [0, v_mag, 0]

# COE
coe_list = coe.rv2coe(earth_mu,r,v,False)
print(coe_list)
[r,v] = coe.coe2rv(earth_mu,coe_list,False)
print(r)
print(v)

# Perturbation
a_j2_6 = perturbations.j2_6(r)
print(a_j2_6)
A = 5.3*2.0
c_d = 2.2
m = 2906
a_drag = perturbations.drag_exp(r,v,A,m,c_d)
print(a_drag)
c_r = 1.2
t = 10000
a_srp = perturbations.srp(t, datetime.datetime.now(), r, A, m, c_r)
print(a_srp)
# Date
jd_now = time_trans.jd(datetime.datetime.now())
print(jd_now)

r_moon = cbe.moon(jd_now)
print(r_moon)
print(np.linalg.norm(r_moon))
a_n_body_moon = perturbations.n_body(t,r,r_moon,4902)
print(a_n_body_moon)
au = 149597870
[r_sun, rtasc, rdecl] = cbe.sun(jd_now)
print(r_sun*au)
print(np.linalg.norm(r_sun*au))
a_n_body_sun = perturbations.n_body(t,r,r_sun*au,132.712e9)
print(a_n_body_sun)

tle = tle_reader.tle_reader('COSMOS 2251 TLE.txt')
print(tle)
sat = satellite.Satellite('COSMOS 2251 TLE.txt')
print(sat.coe)
print(sat.epoch)
print(planet_data.Mercury.radius)
t_span = [0, 30*24*60*60]
y0 = np.concatenate((sat.r,sat.v))
# sol = solve_ivp(eom.two_body, t_span, y0, method='RK45', args=(earth_mu,), rtol=1e-8, atol=1e-8)
# plotter.orbit_plot(earth_radius, sol.y[0:3], sol.t, sat.name)

sol = solve_ivp(eom.cowell, t_span, y0, method='RK45', args=(True, True, True, True, earth_mu, A, m, c_d, c_r, sat.epoch, "Moon",), rtol=1e-8, atol=1e-8)
plotter.orbit_plot(earth_radius, sol.y[0:3], sol.t, sat.name)
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
[r,v] = coe.coe2rv(earth_mu,coe_list[1:7],False)
print(r)
print(v)

# Perturbation
a_j2_6 = perturbations.j2_6(r)
print(a_j2_6, "j2")
A = 5.3*2.0
c_d = 2.2
m = 2906
a_drag = perturbations.drag_exp(r,v,A,m,c_d)
print(a_drag, "drag")
c_r = 1.2
t = 10000
a_srp = perturbations.srp(t, datetime.datetime.now(), r, A, m, c_r)
print(a_srp, "srp")
# Date
jd_now = time_trans.jd(datetime.datetime.now())
print(jd_now)

r_moon = cbe.moon(jd_now)
print(r_moon)
print(np.linalg.norm(r_moon))
a_n_body_moon = perturbations.n_body(r,r_moon,4902)
print(a_n_body_moon, "moon")
au = 149597870
[r_sun, rtasc, rdecl] = cbe.sun(jd_now)
print(r_sun*au)
print(np.linalg.norm(r_sun*au))
a_n_body_sun = perturbations.n_body(r,r_sun*au,132.712e9)
print(a_n_body_sun, "sun")

tle = tle_reader.tle_reader('COSMOS 2251 TLE.txt')
print(tle)
sat = satellite.Satellite('COSMOS 2251 TLE.txt')
print(sat.coe)
print(sat.epoch)
print(planet_data.Mercury.radius)
t_span = [0, 30*24*60*60]
sat = satellite.Satellite('BSAT-3C.txt')
#sat = satellite.Satellite('LEMUR-2 JOEL.txt')
#sat = satellite.Satellite('MOLNIYA 3-50.txt')
print(sat.tle)
print(sat.coe)
print(sat.r)
print(sat.v)
y0 = np.concatenate((sat.r,sat.v))
print(sat.coe[1:7])
t_span = [0, 60*60]
# sol = solve_ivp(eom.two_body, t_span, y0, method='RK45', args=(earth_mu,), rtol=1e-8, atol=1e-8)
# plotter.orbit_plot(earth_radius, sol.y[0:3], sol.t, sat.name)
#sol = solve_ivp(eom.cowell, t_span, y0, method='RK45', args=(earth_mu, True, True, True, True, A, m, c_d, c_r, sat.epoch, ["Moon", "Sun"],), rtol=1e-8, atol=1e-8)
#plotter.coe_plot(sol.y[0:6], sol.t/(24*60*60), "vector", earth_mu, earth_radius, sat.name, deg=True, diff=False)
y0 = sat.coe[1:3] + list(np.deg2rad(sat.coe[3:7]))
sol = solve_ivp(eom.VoP, t_span, y0, method='RK45', args=(earth_mu, True, True, True, True, A, m, c_d, c_r, sat.epoch, ["Moon", "Sun"],), rtol=1e-8, atol=1e-8)
plotter.coe_plot(sol.y[0:6], sol.t/(24*60*60), "coe", earth_mu, earth_radius, sat.name, deg=False, diff=False)
r0 = [20000, -105000, -19000]
v0 = [0.9, -3.4, -1.5]
d_t = 2*60*60
r, v = eom.kepler(r0, v0, d_t, earth_mu)
print(r,v)
d_t = 10
r, v, y, t= eom.encke(t_span, d_t, sat.r, sat.v, earth_mu, earth_radius, 100, 1e-8,
    drag=True, srp=True, n_body=True, j2_6=True,
    **pert_info)
plotter.coe_plot(y, t, "vector", earth_mu, earth_radius, sat.name, deg=True, diff=False)
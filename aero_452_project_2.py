# Project 2
# AERO 452
# Ryo Takatori

import numpy as np
from scipy.integrate import solve_ivp
import time
from orbit_tools import classical_orbital_elements as coe, lambert
from orbit_tools import equation_of_motion as eom
from orbit_tools import perturbations
from orbit_tools import planet_data as p_d
from orbit_tools import plotter
from orbit_tools import utility as util

earth = p_d.Earth()
mue = earth.mu
re = earth.radius
c_d = 2.2
c_r = 1.2
t_span_type = [[0, 30*24*60*60], [0, 365*24*60*60]]
t_span_type = [[0, 30*24*60*60]]

def main():
    sat_dict = {
        "LEMUR-2 JOEL":{
            "name": "LEMUR-2 JOEL.txt",
            "planets": ["Sun", "Moon"],
            "spacecraft": {
                "area": 0.1*0.1*3,
                "mass": 4
            },
            "eom": ["vop"]
        },
        "BSAT-3C":{
            "name": "BSAT-3C.txt",
            "planets": ["Sun", "Moon"],
            "spacecraft": {
                "area": 5.3*2.0,
                "mass": 2906
            },
            "eom": ["vop", "cowell"]
        },
        "MOLNIYA 3-50":{
            "name": "MOLNIYA 3-50.txt",
            "planets": ["Sun", "Moon"],
            "spacecraft": {
                "area": 4.4*1.4*6,
                "mass": 1600
            },
            "eom": ["vop"]
        }
    }
    for key in sat_dict:
        for t_span in t_span_type:
            d = sat_dict[key]
            sat = util.Satellite(d["name"])
            y0 = np.concatenate((sat.r,sat.v))
            perts = perturbations.Perturbation()
            area = d["spacecraft"]["area"]
            mass = d["spacecraft"]["mass"]
            perts.set_all(area, mass, c_d, c_r, sat.epoch, ["Moon", "Sun"], 6)
            unpert = solve_ivp(eom.two_body, [t_span[0], t_span[1]+sat.T/2], y0, method='RK45', args=(mue,), rtol=1e-8, atol=1e-8)
            rf = np.array([comp[-1] for comp in unpert.y[0:3]])
            vf = np.array([comp[-1] for comp in unpert.y[3:6]])
            for eom_method in d["eom"]:
                if eom_method == "vop":
                    tic = time.time()
                    y0 = sat.coe[1:3] + list(np.deg2rad(sat.coe[3:7]))
                    sol = solve_ivp(eom.vop, t_span, y0, method='RK45', events=eom.alt_checker_vop,
                        args=(mue, perts,), rtol=1e-8, atol=1e-8)
                    r0, v0 = coe.coe2rv(mue, [comp[-1] for comp in sol.y[0:6]], deg=False)
                    r0 = np.array(r0)
                    v0 = np.array(v0)
                    toc = time.time()
                    print(f"{sat.name} ({eom_method}) : Elapsed Time {str(toc-tic)} [s]")
                    plotter.coe_plot(sol.y[0:6], sol.t/(24*60*60), "coe",
                        mue, re, sat.name, deg=False, diff=False)
                if eom_method == "cowell":
                    tic = time.time()
                    y0 = np.concatenate((sat.r,sat.v))
                    sol = solve_ivp(eom.cowell, t_span, y0, method='RK45', events=eom.alt_checker,
                        args=(mue, perts,), rtol=1e-8, atol=1e-8)
                    r0 = np.array([comp[-1] for comp in sol.y[0:3]])
                    v0 = np.array([comp[-1] for comp in sol.y[3:6]])
                    toc = time.time()
                    print(f"{sat.name} ({eom_method}) : Elapsed Time {str(toc-tic)} [s]")
                    plotter.coe_plot(sol.y[0:6], sol.t/(24*60*60), "vector",
                        mue, re, sat.name, deg=True, diff=False)
                if eom_method == "encke":
                    tic = time.time()
                    _, _, y, t= eom.encke(t_span, d["d_t"], sat.r, sat.v, mue, 1e-8, perts)
                    print(y)
                    r0 = np.array([comp[-1] for comp in y[0:3]])
                    v0 = np.array([comp[-1] for comp in y[3:6]])
                    toc = time.time()
                    print(f"{sat.name} ({eom_method}) : Elapsed Time {str(toc-tic)} [s]")
                    plotter.coe_plot(y, t, "vector", mue, re, sat.name, deg=True, diff=False)
                v0_corr, vf_corr = lambert.universal_variable(r0, rf, sat.T/2, mue)
                deltav0 = np.linalg.norm(v0_corr - v0)
                deltavf = np.linalg.norm(vf_corr - vf)
                deltav = deltav0 + deltavf
                print("Delta-V: ", deltav)

if __name__ == "__main__":
    main()
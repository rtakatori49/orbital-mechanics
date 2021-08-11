# Project 1
# AERO 452
# Ryo Takatori

import numpy as np
from orbit_tools import rendezvous

# Variables
zero_array = np.array([0, 0, 0])
one_day = 24*60*60

def main():
    # Define rendezvous
    rend_dict = {
        # 100 [km] => 40 [km]
        "100_40": {
            "r0": np.array([0, 100, 0]),
            "rf": np.array([0, 40, 0]),
            "v0": zero_array,
            "vf": zero_array,
            "t": one_day,
            "maneuver": "hop",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "delta_v": 0,
            "r_rel_traj": [],
            "title": "100 [km] => 40 [km] Hop"
        },
        # 20~40 [km] hold (football)
        "40_football": {
            "r0": np.array([0, 40, 0]),
            "t": one_day,
            "maneuver": "football",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "delta_v": 0,
            "r_rel_traj": [],
            "title": "20~40 [km] hold (football)"
        },
        # 40 [km] => 1 [km]
        "40_1": {
            "r0": np.array([0, 40, 0]),
            "rf": np.array([0, 1, 0]),
            "v0": zero_array,
            "vf": zero_array,
            "t": one_day/2,
            "maneuver": "hop",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "delta_v": 0,
            "r_rel_traj": [],
            "title": "40 [km] => 1 [km] Hop"
        },
        # 1 [km] Hold
        "1_hold": {
            "r0": np.array([0, 1, 0]),
            "t": one_day,
            "maneuver": "hold",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "title": "1 [km] Hold"
        },
        # 1 [km] => 300 [m]
        "1_300": {
            "r0": np.array([0, 1, 0]),
            "rf": np.array([0, 0.3, 0]),
            "v0": zero_array,
            "vf": zero_array,
            "t": one_day/4,
            "maneuver": "hop",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "delta_v": 0,
            "r_rel_traj": [],
            "title": "1 [km] => 300 [m] Hop"
        },
        # 300 [m] Hold
        "300_hold": {
            "r0": np.array([0, 0.3, 0]),
            "t": one_day,
            "maneuver": "hold",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "title": "300 [m] Hold"
        },
        # 300 [m] => 20 [m]
        "300_20": {
            "r0": np.array([0, 0.3, 0]),
            "rf": np.array([0.02, 0, 0]),
            "v0": zero_array,
            "vf": zero_array,
            "t": one_day/4,
            "maneuver": "hop",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "delta_v": 0,
            "r_rel_traj": [],
            "title": "300 [m] => 20 [m] Hop"
        },
        # 20 [m] hold (R-Bar)
        "20_hold": {
            "r0": np.array([0.02, 0, 0]),
            "v0": zero_array,
            "t": one_day,
            "maneuver": "r_bar",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "delta_v": 0,
            "r_rel_traj": [],
            "title": "20 [m] hold (R-Bar)"
        },
        # 20 [m] => 0 [m] (V-Bar)
        "20_20": {
            "r0": np.array([0.02, 0, 0]),
            "rf": np.array([0, 0.02, 0]),
            "v0": zero_array,
            "vf": zero_array,
            "t": one_day/2,
            "maneuver": "hop",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "delta_v": 0,
            "r_rel_traj": [],
            "title": "20 [m] => 0 [m] (V-Bar) Hop"
        },
        # 20 [m] => 10 [m]
        "20_10": {
            "r0": np.array([0, 0.02, 0]),
            "rf": np.array([0, 0.01, 0]),
            "v0": zero_array,
            "vf": zero_array,
            "t": one_day/4,
            "maneuver": "hop",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "delta_v": 0,
            "r_rel_traj": [],
            "title": "20 [m] => 10 [m] Hop"
        },
        # 10 [m] Hold
        "10_hold": {
            "r0": np.array([0, 0.01, 0]),
            "t": one_day/2,
            "maneuver": "hold",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "title": "10 [m] Hold"
        },
        # 10 [m] => 5 [m]
        "10_5": {
            "r0": np.array([0, 0.01, 0]),
            "rf": np.array([0, 0.005, 0]),
            "v0": zero_array,
            "vf": zero_array,
            "t": one_day/4,
            "maneuver": "hop",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "delta_v": 0,
            "r_rel_traj": [],
            "title": "10 [m] => 5 [m] Hop"
        },
        # 5 [m] Hold
        "5_hold": {
            "r0": np.array([0, 0.005, 0]),
            "t": one_day/2,
            "maneuver": "hold",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "title": "5 [m] Hold"
        },
        # Final Approach (V-Bar)
        "5_0": {
            "r0": np.array([0, 0.005, 0]),
            "rf": zero_array,
            "t": one_day,
            "maneuver": "v_bar",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "delta_v": 0,
            "r_rel_traj": [],
            "title": "Final Approach (V-Bar)"
        },
        # Final Hold
        "0_hold": {
            "r0": zero_array,
            "t": one_day,
            "maneuver": "hold",
            "r_start": [],
            "v_start": [],
            "r_end": [],
            "v_end": [],
            "title": "Final Hold"
        },
    }
    rendezvous.rendezvous("INTELSAT_18_TLE.txt", rend_dict)

if __name__ == "__main__":
    main()
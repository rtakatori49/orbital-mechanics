# Project 1
# AERO 452
# Ryo Takatori

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from . import equation_of_motion as eom
from . import planet_data as p_d
from . import utility as util

# Earth information
earth = p_d.Earth()
mue = earth.mu
re = earth.radius

# Return end of solution
def return_end(sol, z):
    r_list = sol.y[z[0]:int(sum(z)/2)] # position vector list
    v_list = sol.y[int(sum(z)/2):z[1]] # velocity vector list
    r_end = np.array([r_list[x][-1] for x in range(3)]) # end position
    v_end = np.array([v_list[x][-1] for x in range(3)]) # end velocity
    return r_end, v_end

# Plot trajectory
def plot_traj(r_rel_traj, title):
    plt.plot(r_rel_traj[1], r_rel_traj[0])
    plt.title(title + " Rendezvous Trajectory")
    plt.xlabel("y [km]")
    plt.ylabel("x [km]")
    plt.show()

# Setup and rk45 for LEOM
def eom_leom(t_span, y0, r0, r_bar, v_bar, v_c, title, plot_bool):
    sol = solve_ivp(eom.leom, t_span, y0, method='RK45', args=(mue, r_bar, v_bar, v_c, t_span[0]), rtol=1e-8, atol=1e-8)
    r_rel_list = sol.y[6:9] # trajectory
    r_end, v_end= return_end(sol, [0, 6]) # end vectors
    r_rel_end, v_rel_end = return_end(sol, [6, 12]) # end relative vectors
    # Plot if true
    if plot_bool:
        # Adjust plot for r0
        if r_bar:
            plot_traj(r_rel_list, title) # plot
        else:
            plot_traj(np.transpose(np.transpose(r_rel_list) + r0), title) # plot
    return r_end, v_end, r_rel_list, r_rel_end, v_rel_end

# Hop/move
def move(r0, rf, v0, vf, n, t, r, v, r_bar, v_bar, v_c, title):
    r_rel = rf - r0 # relative position
    v_rel = vf - v0 # relative velocity
    t_span = [0, t] # time span
    delta_v, delta_v_1, _ = eom.two_impluse(r_rel, v_rel, n, t) # two impulse maneuver
    y0 = np.concatenate((r, v, r_rel, delta_v_1)) # initial states
    r_end, v_end, r_rel_list, _, _ = eom_leom(t_span, y0, r0, r_bar, v_bar, v_c, title, True) # LEOM
    return r_end, v_end, delta_v, r_rel_list

# Hold
def hold(t, r, v):
    t_span = [0, t] # time span
    y0 = np.concatenate((r, v)) # initial states
    sol = solve_ivp(eom.two_body, t_span, y0, method='RK45', args=(mue,), rtol=1e-8, atol=1e-8) # two body
    r_end, v_end= return_end(sol, [0, 6]) # end vectors
    return r_end, v_end

def rendezvous(filename, rend_dict):
    # Variables
    legend_list = [] # plot legend
    total_deltav = [] # total delta-v
    total_time = [] # total time
    prev_key = None

    # Initial Orbit
    sat = util.Satellite(filename) # load TLE
    n = sat.n # n
    sat_title = sat.name # satellite name
    t_span = [0, sat.T] # time span
    y0 = np.concatenate((sat.r,sat.v)) # initial states
    sol = solve_ivp(eom.two_body, t_span, y0, method='RK45', args=(mue,), rtol=1e-8, atol=1e-8) # two body
    
    # Plot orbit
    fig = plt.figure(figsize=(9,9), dpi=100)
    ax = fig.add_subplot(111, projection='3d')
    plt.plot(sol.y[0], sol.y[1], sol.y[2], color='C1')
    # Draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)*re
    y = np.sin(u)*np.sin(v)*re
    z = np.cos(v)*re
    ax.plot_wireframe(x, y, z)    # Axis limits
    max_val = np.max(np.abs(np.linalg.norm(sat.r)))
    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])
    plt.title(sat_title + " Orbit")
    plt.xlabel("x [km]")
    plt.ylabel("y [km]")
    ax.set_zlabel("z [km]")
    plt.show()

    # Run through dictionary
    for key in rend_dict:
        d = rend_dict[key]
        # Make end vectors of previous move, start
        if prev_key is not None:
                 d["r_start"] = rend_dict[prev_key]["r_end"] # start position
                 d["v_start"] = rend_dict[prev_key]["v_end"] # start velocity
        else:
            d["r_start"] = sat.r # start position
            d["v_start"] = sat.v # start velocity
        t = d["t"] # time
        total_time.append(t)
        r_start = d["r_start"] # start position
        v_start = d["v_start"] # start velocity
        r0 = d["r0"] # initial relative positon
        maneuver = d["maneuver"] # maneuver
        title = d["title"] # title
        # Hop
        if maneuver == "hop":
            rf = d["rf"] # final relative position
            v0 = d["v0"] # start relative velocity
            vf = d["vf"] # final relative velocity
            r_end, v_end, delta_v, r_rel_traj = move(r0, rf, v0, vf, n, t, r_start, v_start, False, False, 0, title) # actual move
            # Append to list/dict
            d["delta_v"] = delta_v # delta-v
            total_deltav.append(delta_v) # total delta-v
            d["r_rel_traj"] = np.transpose(np.transpose(r_rel_traj) + r0) # trajectory
        # Hold
        if maneuver == "hold":
            r_end, v_end = hold(t, r_start, v_start) # hold
            d["r_rel_traj"] = r0 # location of hold
        # Football
        if maneuver == "football":
            v_football = (n*40)/2 # velocity
            delta_v = v_football*2 # delta-v required
            v_rel = np.array([v_football, 0, 0]) # relative velocity
            t_span = [0, t] # time span
            y0 = np.concatenate((r_start, v_start, r0, v_rel)) # initial states
            r_end, v_end, r_rel_traj, _, _ = eom_leom(t_span, y0, np.array([0, 0, 0]), False, False, 0, title, True) # actual move
            # Append to list/dict
            d["delta_v"] = delta_v # delta-v
            total_deltav.append(delta_v) # total delta-v
            d["r_rel_traj"] = r_rel_traj # trajectory
        # R-bar or V-bar
        if maneuver == "r_bar" or maneuver == "v_bar":
            # R-bar
            if maneuver == "r_bar":
                v0 = d["v0"]
                v_c = 0
                r_bar = True
                v_bar = False
                y0 = np.concatenate((r_start, v_start, r0, v0))
                plot_bool = False
            # V-bar
            if maneuver == "v_bar":
                rf = d["rf"] # final relative position
                r_rel = rf - r0 # relative position
                v_c = r_rel[1]/t
                v0 = np.array([0, v_c, 0])
                r_bar = False
                v_bar = True
                y0 = np.concatenate((r_start, v_start, r_rel, v0))
                plot_bool = True
            t_span = [0, t] # time span
            _, v_end, _, _, v_rel_end  = eom_leom(t_span, y0, r0, False, False, 0, title, False) # move without the v-bar or r-bar to get delta-v
            delta_v = np.linalg.norm(v_rel_end - v0) # delta-v
            r_end, v_end, r_rel_traj, _, _ = eom_leom(t_span, y0, r0, r_bar, v_bar, v_c, title, plot_bool) # actual move
            # Append to list/dict
            d["delta_v"] = delta_v
            total_deltav.append(delta_v)
            # Adjust plot for r0
            if maneuver == "r_bar":
                d["r_rel_traj"] = r_rel_traj
            if maneuver == "v_bar":
                d["r_rel_traj"] = np.transpose(np.transpose(r_rel_traj) + r0)
        # Append end vectors
        d["r_end"] = r_end # end position vector
        d["v_end"] = v_end # end velocity vector
        prev_key = key # store previous key
    
    # Plot entire maneuver
    for key in rend_dict:
        d = rend_dict[key]
        maneuver = d["maneuver"]
        # Show delta-v for each maneuver
        if maneuver == "hold":
            print(d["title"] + " Delta-V: 0")
        else:
            print(d["title"] + " Delta-V: " + str(d["delta_v"]))
        legend_list.append(d["title"]) # legend
        r = d["r_rel_traj"] # trajectory
        # Hold/R-bar hold
        if maneuver == "hold" or maneuver == "r_bar":
            plt.plot(r[1], r[0], marker=".", markersize=10)
        # All moves
        else:
            plt.plot(r[1], r[0])
    plt.title(sat_title + " Rendezvous Trajectory")
    plt.xlabel("y [km]")
    plt.ylabel("x [km]")
    plt.legend(legend_list)
    plt.show()
    # Show total delta-v and time
    print(sat_title + " Total Delta-V: " + str(sum(total_deltav)*1000) + " [m/s]")
    print(sat_title + " Total Time: " + str(sum(total_time)/(24*60*60)) + " [day]")

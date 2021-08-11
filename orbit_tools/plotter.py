# Orbit Plotter
# Ryo Takatori
# 08/05/2020

import numpy as np
import pandas as pd
from . import classical_orbital_elements as coe
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Orbit Plot with Animation
def orbit_plot(planet_radius, r, t, title):
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(9,9), dpi=100)
    ax = fig.add_subplot(111, projection="3d")

    # Trajectory
    ax.plot(r[0], r[1], r[2], color="#ff7f0e", label="Trajectory")
    ax.plot([r[0][0]], [r[1][0]], [r[2][0]], "o", color="w", markersize=12, label="Initial Position")
    ax.plot([r[0][-1]], [r[1][-1]], [r[2][-1]], "o", color="#17becf", markersize=12, label="Final Position")
    orbit_plot, = ax.plot([], [], [], "o", color="#2ca02c", markersize=12, label="Object")
    
    # Draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)*planet_radius
    y = np.sin(u)*np.sin(v)*planet_radius
    z = np.cos(v)*planet_radius
    ax.plot_wireframe(x, y, z, color="#1f77b4")

    # Inertial arrow
    l = planet_radius*1.5
    x,y,z = [[0,0,0],[0,0,0],[0,0,0]]
    u,v,w = [[l,0,0],[0,l,0],[0,0,l]]
    ax.quiver(x,y,z,u,v,w, color="w")
    # Axis limits
    max_val = np.max(np.abs(r))
    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])

    # make the panes transparent
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    # make the grid lines transparent
    ax.xaxis._axinfo["grid"]["color"] = (1,1,1,0)
    ax.yaxis._axinfo["grid"]["color"] = (1,1,1,0)
    ax.zaxis._axinfo["grid"]["color"] = (1,1,1,0)

    # Axis labels
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.set_zlabel("z [km]")

    # Title
    ax.set_title(f"{title} Orbit")
    plt.legend(loc="upper right")
    
    # Slider
    axcolor = "black"
    axorbit = plt.axes([0.1, 0.05, 0.8, 0.025], facecolor=axcolor) # axis location
    sorbit = Slider(axorbit, "Index", 0, len(r[0])-1, 0, valstep=1, facecolor="#1f77b4") # slider
    
    # Time text
    props = dict(boxstyle="round", facecolor="k", alpha=0.5)
    time_text = plt.text(0.05, 0.05, "Current Time:\n0 [s]\n0 [h]\n0 [d]", horizontalalignment="center", verticalalignment="center", transform=ax.transAxes, bbox=props)

    # Update UI if slider has been moved; adjust image frame accordingly
    def update(val):
        orbit_plot.set_data(r[0][val], r[1][val])
        orbit_plot.set_3d_properties(r[2][val])
        time_text.set_text("Current Time:\n"+str(round(t[val],2))+" [s]\n"+\
            str(round(t[val]/(60*60),2))+" [h]\n"+\
            str(round(t[val]/(60*60*24),2))+" [d]")
        fig.canvas.draw_idle()
    sorbit.on_changed(update)
    plt.show()

# COE plots
def coe_plot(data, t, data_type, mu, radius, title, deg, diff):
    def loop(prev_coe):
        rp = coe_temp["a"]*(1-coe_temp["ecc"]) # perigee [km]
        ra = (2*coe_temp["a"]) - rp # apogee [km]
        coe_temp["zp"] = rp - radius # perigee altitude [km]
        coe_temp["za"] = ra - radius # apogee altitude [km]
        # COE difference calculation
        if diff:
            # First difference as 0
            if prev_coe == None:
                for coe_ind in coe_list.values():
                    coe_ind.append(0)
                prev_coe = 1
                coe_init = coe_temp.copy() # store previous coe
            # Difference calculation
            else:
                for coe_ind, coe_temp_ind, coe_init_ind in zip(coe_list.values(), coe_temp.values(), coe_init.values()):
                    coe_ind.append(coe_temp_ind - coe_init_ind)
        # Raw COE
        else:
            for coe_ind, coe_temp_ind in zip(coe_list.values(), coe_temp.values()):
                coe_ind.append(coe_temp_ind)
    # Create dictionary to store
    coe_list = {
        "a": [],
        "zp": [],
        "za": [],
        "ecc": [],
        "h": [],
        "inc": [],
        "raan": [],
        "w": [],
        "theta": []
    }
    coe_temp = {
        "a": 0,
        "zp": 0,
        "za": 0,
        "ecc": 0,
        "h": 0,
        "inc": 0,
        "raan": 0,
        "w": 0,
        "theta": 0
    }
    if data_type == "coe":
        orbit_x = []
        orbit_y = []
        orbit_z = []
        for coe_temp["ecc"], coe_temp["h"], coe_temp["inc"], coe_temp["raan"], coe_temp["w"], coe_temp["theta"] in zip(data[0], data[1], data[2], data[3], data[4], data[5]):
            coe_temp["a"] = (coe_temp["h"]**2)/(mu*(1 - coe_temp["ecc"]**2))
            prev_coe = None
            loop(prev_coe)
            r, _ = coe.coe_to_state(mu, [coe_temp["ecc"], coe_temp["h"], coe_temp["inc"], coe_temp["raan"], coe_temp["w"], coe_temp["theta"]], deg=False)
            r.tolist()
            orbit_x.append(r[0])
            orbit_y.append(r[1])
            orbit_z.append(r[2])
    elif data_type == "vector":
        # Preallocate r and v vector
        r = [0, 0 ,0]
        v = [0, 0 ,0]
        prev_coe = None
        for r[0], r[1], r[2], v[0], v[1], v[2] in zip(data[0], data[1], data[2], data[3], data[4], data[5]):
            # Conver r and v vector to COE
            coe_temp["a"], coe_temp["ecc"], coe_temp["h"], coe_temp["inc"], coe_temp["raan"], coe_temp["w"], coe_temp["theta"] = coe.state_to_coe(mu, r, v, deg=deg)
            loop(prev_coe)
    else:
        print("incorrect type")

    # Plot
    data_list = [coe_list["ecc"], coe_list["inc"], coe_list["raan"], coe_list["w"]]
    title_list = ["Altitude", "Eccentricity", "Inclination", "Right Ascension of Ascending Node", "Argument of Periapse"]
    if deg:
        ylabel_list = ["Altitude [km]", "ecc", "inc [deg]", "$\Omega$ [deg]", "$\omega$ [deg]"]
    else:
        ylabel_list = ["Altitude [km]", "ecc", "inc [rad]", "$\Omega$ [rad]", "$\omega$ [rad]"]
    fig = plt.figure()
    for x in range(5):
        plt.subplot(2,3,x+1)
        # Altitude
        if x == 0:
            plt.plot(t, coe_list["zp"], t, coe_list["za"])
            plt.legend(["Periapse", "Apoapse"])
        # Ecc, inc, raan, w
        else:
            plt.plot(t, data_list[x-1])
        plt.title(title_list[x])
        plt.xlabel("Time [day]")
        plt.ylabel(ylabel_list[x])
        # Orbit
        ax = fig.add_subplot(2,3,6, projection="3d")
    if data_type == "vector":
        plt.plot(data[0], data[1], data[2], color="C1")
        max_val = np.max(np.abs(data[0:2]))
    else:
        plt.plot(orbit_x, orbit_y, orbit_z, color="C1")
        max_val = np.max(np.abs([orbit_x, orbit_y, orbit_z]))
    # Draw sphere
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = np.cos(u)*np.sin(v)*radius
    y = np.sin(u)*np.sin(v)*radius
    z = np.cos(v)*radius
    ax.plot_wireframe(x, y, z)    # Axis limits
    ax.set_xlim([-max_val, max_val])
    ax.set_ylim([-max_val, max_val])
    ax.set_zlim([-max_val, max_val])
    plt.title("Orbit")
    plt.xlabel("x [km]")
    plt.ylabel("y [km]")
    ax.set_zlabel("z [km]")
    plt.suptitle(title + " COE Over Time")
    plt.show()

def coe_table(data):
    coe_data = {
        "Values": data,
        "Units": ["[km]", " ", "[km^2/s]", "[deg]", "[deg]", "[deg]", "[deg]"]
        }
    column_name = ["Values", "Units"]
    index_name = ["Semi-Major Axis (a)", "Eccentricity (ecc)", \
        "Angular Momentum (h)", "Inclination (inc)", \
        "Right Ascension of Ascending Node (raan)", "Argument of Periapse (w)", \
        "True Anomaly (theta)"]
    df = pd.DataFrame(coe_data, columns=column_name, index=index_name)
    print(df)
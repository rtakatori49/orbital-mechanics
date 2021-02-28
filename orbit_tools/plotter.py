# Orbit Plotter
# Ryo Takatori
# 08/05/2020

import numpy as np
import PIL
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# Orbit Plot with Animation
def orbit_plot(planet_radius, r, t, title):
    plt.style.use('dark_background')
    fig = plt.figure(figsize=(9,9), dpi=100)
    ax = fig.add_subplot(111, projection='3d')

    # Trajectory
    ax.plot(r[0], r[1], r[2], color='#ff7f0e', label='Trajectory')
    ax.plot([r[0][0]], [r[1][0]], [r[2][0]], 'o', color='w', markersize=12, label='Initial Position')
    ax.plot([r[0][-1]], [r[1][-1]], [r[2][-1]], 'o', color='#17becf', markersize=12, label='Final Position')
    orbit_plot, = ax.plot([], [], [], 'o', color='#2ca02c', markersize=12, label='Object')
    
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
    ax.quiver(x,y,z,u,v,w, color='w')
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
    ax.xaxis._axinfo["grid"]['color'] = (1,1,1,0)
    ax.yaxis._axinfo["grid"]['color'] = (1,1,1,0)
    ax.zaxis._axinfo["grid"]['color'] = (1,1,1,0)

    # Axis labels
    ax.set_xlabel('x [km]')
    ax.set_ylabel('y [km]')
    ax.set_zlabel('z [km]')

    # Title
    ax.set_title(f'{title} Orbit')
    plt.legend(loc="upper right")
    
    # Slider
    axcolor = 'black'
    axorbit = plt.axes([0.1, 0.05, 0.8, 0.025], facecolor=axcolor) # axis location
    sorbit = Slider(axorbit, 'Index', 0, len(r[0])-1, 0, valstep=1, facecolor="#1f77b4") # slider
    
    # Time text
    props = dict(boxstyle='round', facecolor='k', alpha=0.5)
    time_text = plt.text(0.05, 0.05, "Current Time:\n0 [s]\n0 [h]\n0 [d]", horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, bbox=props)

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
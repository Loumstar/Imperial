import numpy as np

from matplotlib import pyplot as plot
from mpl_toolkits import mplot3d as plot3d

import time
import types

"""
class Bar:
    def __init__(self, length, thermal_diffusivity, numberof_gp):
        self.length = length
        self.a = thermal_diffusivity
        self.numberof_gp = numberof_gp

        self.boundary_conditions = []

    def __set_boundary_conditions(self, fnc, i):
        if type(fnc) != types.FunctionType:
            return TypeError("First argument must be a function")
        self.boundary_conditions.append((i, fnc))

    def set_neumann_boundary

    def

def animate_temperature_distribution(temperature_distribution, figure, axis):
    figure.show()
    figure.canvas.draw()

    line, = axis.plot([], [])
    plot.pause(0.5)
    for temperature_snapshot in temperature_distribution:
        x_values, temp_values = unpack_coords(temperature_snapshot)
        
        line.set_data(x_values, temp_values)

        axis.relim()
        axis.autoscale_view(True, True, True)

        figure.canvas.draw()

        print(temp_values[:5])
"""

def unpack_coords(coords):
    list_of_lists = []
    
    if len(coords) != 0:
        no_of_values_to_unpack = len(coords[0])
        
        for n in range(no_of_values_to_unpack):
            list_of_lists.append(list(map(lambda x: x[n], coords)))
    
    return tuple(list_of_lists)

def next_point_temperature(t0, t1, t2, dx, dt, a):
    first_term = a * dt * (dx ** -2) * (t0 + t2)
    second_term = t1 * (1 - (2 * a * dt * (dx ** -2)))
    
    return first_term + second_term

def get_temperature_distribution(x_limits, time_end, bar_temperature, inital_temperature_distribution, numberof_ti, numberof_gp):
    dx = (x_limits[1] - x_limits[0]) / numberof_gp
    dt = time_end / numberof_ti
    phi = [inital_temperature_distribution]
    
    for p in range(1, numberof_ti):
        row = []
        for i in range(numberof_gp):
            row.append(bar_temperature(phi[p - 1], i, dx, dt, numberof_gp))
        phi.append(row)

    return phi

def is_stable(dt, dx, a):
    return a * dt * (dx ** -2) < 0.5

def next_2d_point_temperature(tx0y1, tx1y0, tx1y1, tx1y2, tx2y1, ds, dt, a):
    first_term = (tx2y1 - (2 * tx1y1) + tx0y1) / (ds ** 2)
    second_term = (tx1y2 - (2 * tx1y1) + tx1y0) / (ds ** 2)

    return tx1y1 + (a * dt * (first_term + second_term))

def get_2d_temperature_distribution(x_limits, y_limits, time_end, surface_temperature, inital_temperature_distribution, numberof_ti, numberof_axis_gp):
    dx = (x_limits[1] - x_limits[0]) / numberof_gp
    dt = time_end / numberof_ti
    phi = [inital_temperature_distribution]
    
    for p in range(1, numberof_ti):
        frame = []
        for j in range(numberof_axis_gp):
            row = []
            
            for i in range(numberof_gp):
                row.append(surface_temperature(phi[p - 1], i, j, dx, dt, numberof_gp))
            
            frame.append(row)
        
        phi.append(frame)

    return phi

# Task A

def bar1(t, i, dx, dt, numberof_gp):
    a = 1.172 * (10 ** -5)
    if i == 0 or i == numberof_gp - 1:
        return 323.15
    else:
        return next_point_temperature(t[i - 1], t[i], t[i + 1], dx, dt, a)

numberof_gp = 50
numberof_ti = 3600    

bar1_t0 = [283.15 for n in range(numberof_gp)] # uniform temperature of 10 deg C

bar_temperature_distribution = get_temperature_distribution((0, 0.5), 3600, bar1, bar1_t0, numberof_ti, numberof_gp)

x_values = np.linspace(0, 0.5, num=numberof_gp, endpoint=True)
time_values = np.linspace(0, 3600, num=numberof_ti, endpoint=True)

x_axis, time_axis = np.meshgrid(x_values, time_values)

task_a = plot.figure()
ax1 = task_a.add_subplot(111, projection='3d')
task_a.canvas.set_window_title("Task A")

ax1.plot_surface(x_axis, time_axis, np.array(bar_temperature_distribution))

# Task B

def boundary1(t):
    k = 40
    h = 500
    dx = 0.01
    t_w = 278.15
    return ((h * t_w) + (k * t[1] / dx)) / (h + (k / dx))

def boundary2(t):
    k = 40
    h = 500
    dx = 0.01
    t_w = 278.15
    return ((h * t_w) + (k * t[-2] / dx)) / (h + (k / dx))

def bar2(t, i, dx, dt, numberof_gp):
    a = 1.172 * (10 ** -5)
    if i == 0:
        return boundary1(t)
    elif i == numberof_gp - 1:
        return boundary2(t)
    elif i == numberof_gp // 2:
        return 373.15
    else:
        return next_point_temperature(t[i - 1], t[i], t[i + 1], dx, dt, a)    

bar2_t0 = [373.15 if n == numberof_gp // 2 else 283.15 for n in range(numberof_gp)]

k = 40
h = 500
dx = 0.01
t_w = 278.15

new_bar_temperature_distribution = get_temperature_distribution((0, 0.4), 3600, bar2, bar2_t0, numberof_ti, numberof_gp)

task_b = plot.figure()

ax2 = task_b.add_subplot(111, projection='3d')
task_b.canvas.set_window_title("Task B")

ax2.plot_surface(x_axis, time_axis, np.array(new_bar_temperature_distribution))

# Task C (cba)

# Task D

def surface1(t, i, j, ds, dt, numberof_axis_gp):
    if i == 0 or j == 0:
        return 453.15
    elif i == numberof_axis_gp - 1 or j == numberof_axis_gp - 1:
        return 453.15
    
    if i > 17 and i < 23 and j > 17 and j < 17:
        a = 1.3 * (10 ** -7)
    else:
        a = 1.9 * (10 ** -5)

    return next_2d_point_temperature(t[j][i - 1], t[j - 1][i], t[j][i], t[j + 1][i], t[j][i + 1], ds, dt, a)
    
def surface_t0(numberof_axis_gp):
    t = []
    for j in range(numberof_axis_gp):
        row = []
        for i in range(numberof_axis_gp):
            if (i >= 17 and i <= 23) and (j >= 17 and j <= 23):
                row.append(258.15)
            else:
                row.append(298.15)
        t.append(row)
    return t

k = 40
h = 500
dx = 0.01
t_w = 278.15
numberof_axis_gp = 40

x_values = np.linspace(0, 0.4, num=numberof_gp, endpoint=True)
y_values = np.linspace(0, 0.4, num=numberof_gp, endpoint=True)

x_axis, y_axis = np.meshgrid(x_values, y_values)

surface_temperature_distribution = get_2d_temperature_distribution((0, 0.4), (0, 0.4), 3600, surface1, surface_t0(numberof_axis_gp), numberof_ti, numberof_axis_gp)

task_d1 = plot.figure()

ax3 = task_d1.add_subplot(111, projection='3d')
task_d1.canvas.set_window_title("Task D Part A")

ax3.plot_surface(x_axis, y_axis, np.array(surface_temperature_distribution[0]))

task_d2 = plot.figure()

print(surface_temperature_distribution[3599])

ax4 = task_d2.add_subplot(111, projection='3d')
task_d2.canvas.set_window_title("Task D Part B")

ax4.plot_surface(x_axis, y_axis, np.array(surface_temperature_distribution[3599]))

plot.show()
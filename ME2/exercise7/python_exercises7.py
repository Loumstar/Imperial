import math as m

from matplotlib import pyplot as plot

def unpack_coords(coords):
    list_of_lists = []
    
    if len(coords) != 0:
        no_of_values_to_unpack = len(coords[0])
        
        for n in range(no_of_values_to_unpack):
            list_of_lists.append(list(map(lambda x: x[n], coords)))
    
    return tuple(list_of_lists)

def plot_coords(coords, axis):
    unpacked_coords = unpack_coords(coords)
    domain_values = unpacked_coords[0]
    
    for axis_values in unpacked_coords[1:]:
        axis.plot(domain_values, axis_values)

def distribution(fnc, x_lb, x_ub, h):
    coords = []
    x = x_lb
    
    while x <= x_ub:
        coords.append(tuple([x, fnc(x)]))
        x += h
    
    return coords

def forward_euler_approximation(fnc, lb_coord, x_ub, h):
    coords = []
    x, y = lb_coord
    
    while x < x_ub:
        y = y + (h * fnc(x, y))
        coords.append(tuple([x, y]))
        x += h
    
    return coords

def vector_forward_euler_approximation(fnc, lb_coord, x_ub, h):
    coords = []
    coord = list(lb_coord)
    
    while coord[0] < x_ub:
        vector_derivative = fnc(*coord)
        
        for n in range(len(vector_derivative)):
            coord[n] += h if n == 0 else h * vector_derivative[n]   
        
        coords.append(tuple(coord))

    return coords

def runge_kutta_k_values(n, fnc, coord, h):
    x, y = coord
    k_values = []
    k = 0
    for i in range(n):
        if i == 0:
            k = h * fnc(x, y)
        elif i == n - 1:
            k = h * fnc(x + h, y + k)
        else:
            k = h * fnc(x + (h / 2), y + (k / 2))
        
        k_values.append(k)

    return k_values

def adaptive_simpson_approximation(y_values, h):
    summation = 0
    for i, y in enumerate(y_values):
        if (i == 0) or (i == len(y_values) - 1):
            summation += y
        elif i % 2 == 0:
            summation += 2 * y
        else:
            summation += 4 * y

    return summation * h / 3

def runge_kutta_approximation(n, fnc, lb_coord, x_ub, h):
    coords = []
    x, y = lb_coord
    
    while x < x_ub:
        k_values = runge_kutta_k_values(n, fnc, (x, y), h)
        y += adaptive_simpson_approximation(k_values, h) / (2 * h)
        coords.append(tuple([x, y]))
        x += h
    
    return coords

def backward_euler_approximation(fnc, x_lb, ub_coord, h):
    coords = []
    x, y = ub_coord
    
    while x > x_lb:
        y = y - (h * fnc(x, y))
        coords.append(tuple([x, y]))
        x -= h
    
    return coords

"""
Main Script
"""

# Task A

def f1(t, y):
    return -2 * ((y * t) + (t ** 3))

def f2(t):
    return 1 - (t ** 2)

f1_analytical_solution = distribution(f2, 0, 5.0, 0.01)

f1_forward_euler_integral = forward_euler_approximation(f1, (0, 1), 5.0, 0.1)
f1_runge_kutta_integral = runge_kutta_approximation(4, f1, (0, 1), 5.0, 0.1)

# Task B

#f1_backward_euler_integral = backward_euler_approximation(f1, 0, (5.0, -24.3), 0.1)

task_ab, ax1 = plot.subplots(1, 1)
task_ab.canvas.set_window_title("Task A and B")

plot_coords(f1_analytical_solution, ax1)
plot_coords(f1_forward_euler_integral, ax1)
plot_coords(f1_runge_kutta_integral, ax1)
#plot_coords(f1_backward_euler_integral, ax1)

# Task C

def rent_price_rate(bps, p):
    return bps * ((0.3 * p) - 0.8)

def population_rate(bps, p):
    return p * (1.1 - bps)

def rental_rates(t, bps, p):
    return 0, rent_price_rate(bps, p), population_rate(bps, p)

rental_variation = vector_forward_euler_approximation(rental_rates, (0, 0.8, 7), 40, 0.019)

task_c, (ax2, ax3) = plot.subplots(2, 1)
task_c.canvas.set_window_title("Task C")
plot_coords(rental_variation, ax2)

rent_price_data, population_data = unpack_coords(rental_variation)[1:]
ax3.plot(rent_price_data, population_data)

# Task D

def spring_velocity(t, theta, w):
    return w

def dry_spring_acceleration(t, theta, w):
    damping_coeff = 0.05
    mass = 0.50
    length = 1
    g = 9.81
    return -((damping_coeff * w / mass) + (g * m.sin(theta) / length))

def wet_spring_acceleration(t, theta, w):
    damping_coeff = 0.18
    mass = 0.50
    length = 1
    g = 9.81
    return -((damping_coeff * w / mass) + (g * m.sin(theta) / length))

def dry_spring_rates(t, theta, w):
    return 0, spring_velocity(t, theta, w), dry_spring_acceleration(t, theta, w)

def wet_spring_rates(t, theta, w):
    return 0, spring_velocity(t, theta, w), wet_spring_acceleration(t, theta, w)

dry_pendulum_motion = vector_forward_euler_approximation(dry_spring_rates, (0, m.pi / 4, 0), 15, 0.005)
wet_pendulum_motion = vector_forward_euler_approximation(wet_spring_rates, (0, m.pi / 4, 0), 15, 0.005)

task_d, (ax4, ax5) = plot.subplots(2, 1)
task_d.canvas.set_window_title("Task D")

plot_coords(dry_pendulum_motion, ax4)
plot_coords(wet_pendulum_motion, ax5)

plot.show()
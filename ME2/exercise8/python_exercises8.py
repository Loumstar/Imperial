import numpy as np
import math as m

from matplotlib import pyplot as plot

def solve_ode_direct(fnc, x_limits, b_values, b_coeff, numberof_gp):
    x_lb, x_ub = x_limits
    b_lb, b_ub = b_values

    matrix_of_constants = []

    x_values = []
    p_values = []

    step = (x_ub - x_lb) / (numberof_gp)

    for n in range(numberof_gp):
        row = np.zeros(numberof_gp)
        x = x_lb + (n * step)

        if n == 0:
            p_values.append(b_lb)
            row[0] = b_coeff[1] - (b_coeff[0] / step)
            row[1] = b_coeff[0] / step
        elif n == numberof_gp - 1:
            p_values.append(b_ub)
            row[n] = 1
            row[n - 1] = -b_coeff[2] / step
            row[n] = (b_coeff[2] / step) + b_coeff[3]
        else:
            f, g, h = fnc(x)
            row[n - 1] = (step ** -2) - (f / (2 * step))
            row[n] = g - (2 * (step ** -2))
            row[n + 1] = (step ** -2) + (f / (2 * step))
            
            p_values.append(h)

        matrix_of_constants.append(row)
        x_values.append(x)

    y_values = np.linalg.inv(matrix_of_constants).dot(p_values)

    return x_values, y_values

def jacobi_approximation(y_n_minus_1, y_n_plus_1, f, g, h, step):
    y_n_plus_1_term = y_n_plus_1 * ((step ** -2) + (f / (2 * step)))
    y_n_minus_1_term = y_n_minus_1 * ((step ** -2) - (f / (2 * step)))
    y_n_coefficient = g - (2 * (step ** -2))
    
    return (h - y_n_plus_1_term - y_n_minus_1_term) / y_n_coefficient

def jacobi_next_iteration(fnc, x, y, numberof_gp, step):
    new_y = y.copy()
    
    for n in range(1, numberof_gp - 1):
        f, g, h = fnc(x[n])
        new_y[n] = jacobi_approximation(new_y[n - 1], new_y[n + 1], f, g, h, step)
    
    return new_y


def max_difference(x_values, y_values):
    return max(map(lambda x, y: abs(x - y), x_values, y_values))

def solve_ode_jacobi(fnc, x_limits, b_values, numberof_gp, max_iterations, tolerance):
    x_values = np.linspace(x_limits[0], x_limits[1], numberof_gp, endpoint=True)
    step = (x_limits[1] - x_limits[0]) / (numberof_gp)

    previous_y_values = np.zeros(numberof_gp)
    y_values = np.zeros(numberof_gp)

    y_values[0], y_values[numberof_gp - 1] = b_values

    for _ in range(max_iterations):
        y_values = jacobi_next_iteration(fnc, x_values, y_values, numberof_gp, step)
        if max_difference(y_values, previous_y_values) < tolerance:
            break
        else:
            previous_y_values = y_values.copy()

    return x_values, y_values

# Task A

def ode1(x):
    """
    Returns f(x), g(x), h(x) for the differential equation
    y''(x) + 2x y'(x) + 2y - cos(3x) = 0
    Where:
    y''(x) + f(x) y'(x) + g(x) y(x) = h(x)
    """
    return (2 * x), 2, m.cos(3 * x)

task_ab, ax1 = plot.subplots(1, 1)
task_ab.canvas.set_window_title("Task A and B")

ax1.plot(*solve_ode_direct(ode1, (0, m.pi), (1.5, 0), (0, 1, 0, 1), 100))

# Task B

ax1.plot(*solve_ode_jacobi(ode1, (0, m.pi), (1.5, 0), 100, 1000, 0.0001))

# Task C (cba)

def ode2(x):
    return x, 1, -5 * x

def ode3(x):
    return x, 1, 0 * x

def ode4(x):
    return x, 1, 5 * x

task_c, ax2 = plot.subplots(1, 1)
task_c.canvas.set_window_title("Task C")

ax2.plot(*solve_ode_direct(ode2, (0, 2.0), (0.0, -1.0), (1, 0, 1, 0), 50))
ax2.plot(*solve_ode_direct(ode3, (0, 2.0), (0.0, -1.0), (1, 0, 1, 0), 50))
ax2.plot(*solve_ode_direct(ode4, (0, 2.0), (0.0, -1.0), (1, 0, 1, 0), 50))

# Task D

task_d, ax3 = plot.subplots(1, 1)
task_d.canvas.set_window_title("Task D")

def temperature_ode(r):
    """
    Returns f(x), g(x), h(x) for the differential equation
    T''(r) + (1/r) T'(r) = -10^8 * e^(-r/R) / kr
    Where:
    y''(x) + f(x) y'(x) + g(x) y(x) = h(x)
    """
    R = 0.015
    k = 16.75
    return 1 / r, 0, -(10 ** 8) * (m.e ** (-r / R)) / (k * r)

radius = 0.015
w = 0.003
h = 6 * (10 ** 4)
k = 16.75
t_w = 473

x_coords = (radius, radius + w)
b_values = (-6.32 * (10 ** 5) / k, h * t_w / k)

boundary_coefficients = (1, 0, 1, h / k)

ax3.plot(*solve_ode_direct(temperature_ode, x_coords, b_values, boundary_coefficients, 50))

plot.show()
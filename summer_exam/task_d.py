import math
from matplotlib import pyplot as plot

def integrate_numerically(fnc, x_lb, x_ub, h):
    integral_sum = 0
    x = x_lb
    
    while x + h < x_ub:
        integral_sum += h * (fnc(x) + fnc(x + h)) / 2
        x += h
    
    return integral_sum

def green_region_function(x, n):
    """
    Method that descibes the vertical height in the green region as a function of x
    """
    circle_height  = math.sqrt((10 ** 2) - ((x - 3) ** 2))
    red_function = 5 * math.sin(2 * math.pi * n * x / 13) * (math.e ** -(x / 10))

    return circle_height - red_function

"""
Main Script
"""

cid_values = [0, 1, 3, 6, 0, 4, 5, 8]
integral_values = []

for n in cid_values:
    # Redefine function for different values of n.
    def fnc(x):
        return green_region_function(x, n)

    # Run integral and add to list
    area = integrate_numerically(fnc, 0, 13, 0.02)
    integral_values.append(area)

# Plot list again the values of n
plot.scatter(cid_values, integral_values)
plot.show()
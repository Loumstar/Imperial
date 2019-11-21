import numpy
from matplotlib import pyplot as plot

# Task A

def integrate_numerically(coords):
    x0, y0 = coords[0]
    integral_sum = 0
    
    for x1, y1 in coords[1:]:
        width = x1 - x0
        
        if width < 0:
            print("Coords must be in increasing order.")
            print("%f > %f" % (x0, x1))
            return None

        integral_sum += ((y0 + y1) / 2) * width
        x0, y0 = x1, y1
    
    return integral_sum

def get_function_distribution(fnc, x_lb, x_ub, step):
    distribution = []
    x = x_lb
    while x < x_ub:
        distribution.append((x, fnc(x)))
        x += step
    return distribution

def split_coords(coords):
    x_range, y_range = [], []
    for x, y in coords:
        x_range.append(x)
        y_range.append(y)
    return x_range, y_range

def f1(x):
    return ((x ** 8.11) + 2019) ** -0.5

f1_distribution = get_function_distribution(f1, 0, 2, 0.4)
f1_integral_approximation = integrate_numerically(f1_distribution)

x_range, y_range = split_coords(f1_distribution)

# Task B

def f2(x):
    return ((x ** 1.11) + 2019) ** -0.5

f2_distribution = get_function_distribution(f2, 0, 2, 0.4)
f2_integral_approximation = integrate_numerically(f2_distribution)

# Task C

"""
Already done as part of Task A
"""

# Task D

def join_coords(x_coords, y_coords):
    coords = []

    if len(x_coords) != len(y_coords):
        print("Number of x coords not the same as number of y coords")
    
    numberof_coords = len(x_coords) \
        if len(x_coords) == len(y_coords) \
        or len(x_coords) < len(y_coords) \
        else len(y_coords)
    
    for i in range(numberof_coords):
        coords.append((x_coords[i], y_coords[i]))

    return coords

def integrate_numerically_allow_overlap(coords):
    """
    Method that doesn't check that x is always increasing.
    
    This is because the thames goes back on itself so integrating causes this error in 
    the other method. This will compute the correct area as long as the coords are in 
    order of the boundary path.
    """
    x0, y0 = coords[0]
    integral_sum = 0
    
    for x1, y1 in coords[1:]:
        integral_sum += ((y0 + y1) / 2) * (x1 - x0)

        x0, y0 = x1, y1
    
    return integral_sum

with open('./xn.txt', 'r') as f:
    thames_north_edge_x = [float(xn) for xn in f.read().split('\n')]

with open('./yn.txt', 'r') as f:
    thames_north_edge_y = [float(yn) for yn in f.read().split('\n')]

with open('./xs.txt', 'r') as f:
    thames_south_edge_x = [float(xs) for xs in f.read().split('\n')]

with open('./ys.txt', 'r') as f:
    thames_south_edge_y = [float(ys) for ys in f.read().split('\n')]

plot.plot(
    thames_north_edge_x, thames_north_edge_y, 'r',
    thames_south_edge_x, thames_south_edge_y, 'b'
)

plot.axis('equal')
plot.show()

thames_north_edge_coords = join_coords(thames_north_edge_x, thames_north_edge_y)
thames_south_edge_coords = join_coords(thames_south_edge_x, thames_south_edge_y)

thames_area = integrate_numerically_allow_overlap(thames_north_edge_coords) \
            - integrate_numerically_allow_overlap(thames_south_edge_coords)

print("Area of thames: %.f m2" % thames_area)

# Task E

"""
CBA
"""
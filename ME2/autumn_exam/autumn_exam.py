import math, numpy
from matplotlib import pyplot as plot

# Question 3

def upper_riemann_integral(coords):
    x0, y0 = coords[0]
    integral_sum = 0
    
    for x1, y1 in coords[1:]:
        width = x1 - x0
        
        if width < 0:
            print("Coords must be in increasing order.")
            print("%f > %f" % (x0, x1))
            return None

        upper_value = y1 if y1 > y0 else y0

        integral_sum += upper_value * width
        x0, y0 = x1, y1
    
    return integral_sum

def get_function_distribution(fnc, x_lb, x_ub, step):
    distribution = []
    x = x_lb
    while x < x_ub:
        distribution.append((x, fnc(x)))
        x += step
    return distribution

def f1(x):
    return math.sin(2 * x)

sin2x_coords = get_function_distribution(f1, 6, 11, 0.1)
upper_riemann_sin2x = upper_riemann_integral(sin2x_coords)

print(upper_riemann_sin2x)

# Question 4

def count_pixels_with_red_as(value, image_array):
    count = 0
    for j in range(image_array.shape[0]):
        for i in range(image_array.shape[1]):
            if image_array[j][i][0] == value:
                count += 1
    return count

photo_filename = "NoCheating(3).png"
photo_array = plot.imread(photo_filename).copy()

print(count_pixels_with_red_as(1, photo_array))

# Question 5

def numpy_get_distribution(fnc, x_lb, x_ub, step):
    x_range = numpy.arange(x_lb, x_ub, step)
    coords = [tuple([x, fnc(x)]) for x in x_range]
    return coords

def trapezium_integral(coords):
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

def f2(x):
    return math.e ** -((x ** 2) / 0.8)

f2_coords_1 = numpy_get_distribution(f2, -10, -0.60, 0.05)
f2_coords_2 = numpy_get_distribution(f2, -0.60, 0.61, 0.01)
f2_coords_3 = numpy_get_distribution(f2, 0.65, 10.5, 0.05)

f2_integral_1 = trapezium_integral(f2_coords_1)
f2_integral_2 = trapezium_integral(f2_coords_2)
f2_integral_3 = trapezium_integral(f2_coords_3)

total_f2_integral = f2_integral_1 + f2_integral_2 + f2_integral_3

print(total_f2_integral)

# Question 6

# Run code from exam:
x = numpy.arange(0, 4 * numpy.pi, 0.1)
y = numpy.sin(x)

# x becomes all list of f(x) where the all values of f^2(x) < 2
x = y[ y * y < 2 ]

print(x)

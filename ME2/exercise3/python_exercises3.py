import math, numpy

from matplotlib import pyplot as plot
from mpl_toolkits import mplot3d as plot3d

# From exercise 2

def vector_gradient(vector, x, y, z, wrt, h):
    derivative = [0, 0, 0]
    
    if h:
        p0 = vector(x, y, z)
        
        if wrt == 'x':
            p1 = vector(x + h, y, z)
        elif wrt == 'y':
            p1 = vector(x, y + h, z)
        elif wrt == 'z':
            p1 = vector(x, y, z + h)
        else:
            return (0, 0, 0)

        for d in (0, 1, 2):
            derivative[d] = (p1[d] - p0[d]) / h

    return tuple(derivative)

def curl(vector, x, y, z, h):
    i = vector_gradient(vector, x, y, z, 'y', h)[2] - vector_gradient(vector, x, y, z, 'z', h)[1]
    j = vector_gradient(vector, x, y, z, 'z', h)[0] - vector_gradient(vector, x, y, z, 'x', h)[2]
    k = vector_gradient(vector, x, y, z, 'x', h)[1] - vector_gradient(vector, x, y, z, 'y', h)[0]

    return (round(i, 5), round(j, 5), round(k, 5))

def get_curl_distribution(vector, x_range, y_range, z_range, h):
    vector_field_array = []
    for z in z_range:
        row = []
        for y in y_range:
            column = []
            for x in x_range:
                column.append(curl(vector, x, y, z, h))
            row.append(column)
        vector_field_array.append(row)

    return vector_field_array

# Task A

def scalar_gradient(fnc, x, y, z, wrt, h):
    p0 = fnc(x, y, z)

    if wrt == 'x':
        p1 = fnc(x + h, y, z)
    elif wrt == 'y':
        p1 = fnc(x, y + h, z)
    elif wrt == 'z':
        p1 = fnc(x, y, z + h)
    else:
        return 0

    return (p1 - p0) / h

def sum_differentials(fnc, x_range, y_range, z_range, h):
    sum_array = []
    for z in z_range:
        row = []
        for y in y_range:
            column = []
            for x in x_range:
                diff_sum = \
                    scalar_gradient(fnc, x, y, z, 'x', h) + \
                    scalar_gradient(fnc, x, y, z, 'y', h) + \
                    scalar_gradient(fnc, x, y, z, 'z', h)
                
                column.append(diff_sum)
            row.append(column)
        sum_array.append(row)

    return sum_array

x_range = y_range = z_range = numpy.arange(-2, 2, 0.1)

h = 0.01

def f1(x, y, z):
    return (3 * (x ** 2)) + (math.e ** y) + (z * y)

def v1(x, y, z):
    return 6 * x, z + (math.e ** y), y

f1_diff_sum = sum_differentials(f1, x_range, y_range, z_range, h)
v1_curl_dist = get_curl_distribution(v1, x_range, y_range, z_range, h)

# Task B

def laplacian(fnc, x, y, z, h1, h2):
    # Where h1 << h2
    f_x0 = scalar_gradient(fnc, x, y, z, 'x', h1)
    f_y0 = scalar_gradient(fnc, x, y, z, 'y', h1)
    f_z0 = scalar_gradient(fnc, x, y, z, 'z', h1)

    f_x1 = scalar_gradient(fnc, x + h2, y, z, 'x', h1)
    f_y1 = scalar_gradient(fnc, x, y + h2, z, 'y', h1)
    f_z1 = scalar_gradient(fnc, x, y, z + h2, 'z', h1)

    return ((f_x1 - f_x0) + (f_y1 - f_y0) + (f_z1 - f_z0)) / h2

def f2(x, y, z):
    return (x ** 2) - (y ** 2)

def satisfies_laplace(fnc, x_range, y_range, z_range, h1, h2):
    for z in z_range:
        for y in y_range:
            for x in x_range:
                if round(laplacian(fnc, x, y, z, h1, h2), 5) != 0:
                    return False
                
    return True

does_f2_satisfy_laplace = satisfies_laplace(f2, x_range, y_range, [0], 0.001, 0.1)

print("f2 satisfies the Laplacian: %s." % (does_f2_satisfy_laplace))

# Task C

def get_temperature_distribution(fnc, x_range, y_range, z_range, numberof_terms):
    temperature_distribution_array = []
    for z in z_range:
        row = []
        for y in y_range:
            column = []
            for x in x_range:
                t = 0
                for n in range(numberof_terms):
                    t += fnc(x, y, z, n)
                column.append(t)
            row.append(column)
        temperature_distribution_array.append(row)

    return temperature_distribution_array

def t1(x, y, z, n):
    m = (2 * n) + 1 # N must be odd

    coefficient = 40 / (math.pi * m * math.sinh(-m * math.pi / 2))

    x_function = math.sinh(m * (x - 2) * math.pi / 4)
    y_function = math.sin(m * math.pi * y / 4)

    return coefficient * x_function * y_function

x_range = numpy.arange(0, 2, 0.1)
y_range = numpy.arange(0, 4, 0.1)
z_range = tuple([0])

N = 105

x_axis, y_axis = numpy.meshgrid(x_range, y_range)

caretesian_temp_array = numpy.array(get_temperature_distribution(t1, x_range, y_range, z_range, N)[0])

fig = plot.figure()

temperature_plot = fig.add_subplot(111, projection='3d')
temperature_plot.plot_surface(x_axis, y_axis, caretesian_temp_array)
temperature_plot.axis((0, 2, 0, 4))

# Task D

def to_polar(x, y):
    r = math.hypot(x, y)
    theta = math.atan2(y, x)

    return r, theta

# T = 1 and a = 1
def t2(r, theta, z, n):
    m = (2 * n) + 1

    coefficient = 4 / (math.pi * m)

    r_function = r ** m
    theta_function = math.sin(m * theta)

    return coefficient * r_function * theta_function if r <= 1 else 0

def cartesian_t2(x, y, z, n):
    r, theta = to_polar(x, y)
    return t2(r, theta, z, n)

x_range = numpy.arange(-1, 1, 0.01)
y_range = numpy.arange(0, 1, 0.01)
z_range = tuple([0])

N = 45

x_axis, y_axis = numpy.meshgrid(x_range, y_range)

polar_temp_array = numpy.array(get_temperature_distribution(cartesian_t2, x_range, y_range, z_range, N)[0])

fig = plot.figure()

temperature_plot = fig.add_subplot(111, projection='3d')
temperature_plot.plot_surface(x_axis, y_axis, polar_temp_array)
temperature_plot.axis((-1, 1, 0, 1))

plot.show()

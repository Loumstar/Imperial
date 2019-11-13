from matplotlib import pyplot as plot
from mpl_toolkits import mplot3d as plot3d
import numpy

# Task A

def remove_green(image_array):
    for j in range(image_array.shape[0]):
        for i in range(image_array.shape[1]):
            image_array[j][i][1] = 0
    
    return image_array

def mirror(image_array):
    for j in range(image_array.shape[0]):
        image_array[j] = image_array[j][::-1]
    
    return image_array

def shrink(image_array, n):
    small_array = []
    
    for j in range(int(image_array.shape[0] / n)):
        row = []
        for i in range(int(image_array.shape[1] / n)):
            row.append(image_array[int(n * j)][int(n * i)])
        small_array.append(row)

    return numpy.array(small_array)

parrots_filename = "Parrots.jpg"

parrots_image = plot.imread(parrots_filename).copy()

parrots_image_no_green = remove_green(parrots_image)
parrots_image_no_green_mirrored = mirror(parrots_image_no_green)
parrots_image_no_green_mirrored_small = shrink(parrots_image_no_green_mirrored, 2.5)

plot.imsave('parrots_copy.jpg', parrots_image_no_green_mirrored_small)

# Task B

def vector_gradient(vector, x, y, z, wrt, h):
    derivative = [0, 0, 0]
    
    if h:
        if wrt == 'x':
            p1 = vector(x + h, y, z)
        elif wrt == 'y':
            p1 = vector(x, y + h, z)
        elif wrt == 'z':
            p1 = vector(x, y, z + h)

        p0 = vector(x, y, z)

        for d in (0, 1, 2):
            derivative[d] = (p1[d] - p0[d]) / h

    return tuple(derivative)


def div(vector, x, y, z, h):
    i = vector_gradient(vector, x, y, z, 'x', h)[0]
    j = vector_gradient(vector, x, y, z, 'y', h)[1]
    k = vector_gradient(vector, x, y, z, 'z', h)[2]
    
    print(i, j, k)

    return round(i + j + k, 5)

def get_div_distribution(vector, x_range, y_range, z_range, h):
    scalar_field_array = []
    for z in z_range:
        row = []
        for y in y_range:
            column = []
            for x in x_range:
                column.append(div(vector, x, y, z, h))
            row.append(column)
        scalar_field_array.append(row)

    return scalar_field_array

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

def vector1(x, y, z):
    return (x, y, 0)

def vector2(x, y, z):
    return (y, -x, 0)

def vector3(x, y, z):
    i = y / ((x ** 2) + (y ** 2))
    j = -x / ((x ** 2) + (y ** 2))
    k = 0
    return (i, j, k)

x_range = numpy.arange(-5, 5, 0.1)
y_range = numpy.arange(-5, 5, 0.1)
z_range = [0]

#print(div(vector3, 0.1, 0.2, 0, 0.0001))

"""
div_values = numpy.array(get_div_distribution(vector3, x_range, y_range, z_range, 0.001)[0])
curl_values = numpy.array(get_curl_distribution(vector3, x_range, y_range, z_range, 0.001)[0])

fig = plot.figure()

div_plot = fig.add_subplot(111, projection='3d')
x_axis, y_axis = numpy.meshgrid(x_range, y_range)

div_plot.plot_surface(x_axis, y_axis, div_values)

plot.show()


No code to plot curl vector field
"""

# Task C

def streamline(velocity_vector_2d, x_range, initial_position):
    y = 0
    y_range = [0]
    initial_position_index = None

    for i, x in enumerate(x_range[1:]):
        vector = velocity_vector_2d(x, y)
        dy = (vector[1] / vector[0]) * (x - x_range[i])
        y += dy
        y_range.append(y)

        if initial_position[0] == round(x, 5):
            initial_position_index = i + 1

    if initial_position_index:
        integrand_constant = initial_position[1] - y_range[initial_position_index]
        y_range = list(map(lambda y: y + integrand_constant, y_range))
    else:
        print("Initial position not in the x_range. Integrand constant not found.")

    return y_range

def fluid_velocity(x, y):
    i = (4 * x) + (14 * y)
    j = -(6 * x) - (11 * y)
    return i, j

x_range = numpy.arange(-10, 10, 0.1)
y_range = streamline(fluid_velocity, x_range, (0, 0))

plot.clf()
plot.plot(x_range, y_range, 'r')

plot.show()

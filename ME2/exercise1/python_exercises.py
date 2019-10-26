from matplotlib import pyplot as plot
from math import sin

# Task A

def linear_sequence(lb, ub, step):
    arr = []
    x = lb
    while x < ub:
        arr.append(x)
        x += step
    return arr

x_values = linear_sequence(-5, -2, 0.1)
x_values += linear_sequence(-2, 3, 0.01)
x_values += linear_sequence(3, 5, 0.1)

y_values = [sin(x) for x in x_values]

plot.plot(x_values, y_values, 'r')

plot.show()
plot.clf()

# Task B

def backwards_difference_derivative(x_values, y_values):
    dy_values = []
    
    y0 = y_values[0]
    x0 = x_values[0]

    for i, y1 in enumerate(y_values[1:]):
        x1 = x_values[i + 1]
        dy_values.append((y1 - y0) / (x1 - x0))

        x0, y0 = x1, y1

    return dy_values

def central_difference_derivative(x_values, y_values):
    dy_values = []
    
    y0 = y_values[0]
    x0 = x_values[0]

    for i, y2 in enumerate(y_values[2:]):
        x2 = x_values[i + 2]

        dy_values.append((y2 - y0) / (x2 - x0))

        y0 = y_values[i + 1]
        x0 = x_values[i + 1]

    return dy_values

dy_values_backwards = backwards_difference_derivative(x_values, y_values)

plot.plot(
    x_values, y_values, 'r',
    x_values[1:], dy_values_backwards, 'b'
)

plot.show()
plot.clf()

dy_values_central = central_difference_derivative(x_values, y_values)

plot.plot(
    x_values, y_values, 'r',
    x_values[2:], dy_values_central, 'b'
)

plot.show()
plot.clf()

# Task C

absolute_y_values = [abs(y) for y in y_values]

def y_function(x, y):
    if x <= 0:
        return 0
    elif y <= 0.5:
        return y
    else:
        return 0.5
    
new_y_values = list(map(y_function, x_values, absolute_y_values))

plot.plot(
    x_values, y_values, 'r',
    x_values, new_y_values, 'g'
)

plot.show()

# Task D

def create_3d_array(x, y, z):
    z_arr = []
    for k in range(z):
        y_arr = []
        for j in range(y):
            x_arr = []
            for i in range(x):
                x_arr.append(int((j + k) % 2 == 0))
            y_arr.append(x_arr)
        z_arr.append(y_arr)

    return z_arr

cube = create_3d_array(5, 4, 3)

print(cube)
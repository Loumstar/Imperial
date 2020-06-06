import math, numpy
from matplotlib import pyplot as plot

def vector_forward_euler(derivative_fnc, lb_coord, lb_derivative, x_ub, h):
    """
    Method to return a forward euler approximation for a vector. 
    
    As the first element in the coordinates is the domain variable (x, t, etc.) the first
    element in the tuple returned by the function is the derivative of the domain 
    variable with respect to itself and should therefore be equal to 1.

    lb_coord is the lower bound coordinate vector, and lb_derivative is the lower bound
    derivative vector.
    """
    # Create a list of lists for each set of coordinates
    coords = [[lb_coord[0]], [lb_coord[1]], [lb_coord[2]]]
    i = 0

    while coords[0][i] < x_ub:
        # Calculate the vector derivative using the function that defines the ode
        vector_derivative = derivative_fnc(coords[0][i], coords[1][i], coords[2][i]) \
            if i != 0 else lb_derivative

        # For each new coordinate, apply forwards euler using the vector derivatives
        # and then add the new coordinates to the list
        for n in range(len(vector_derivative)):
            coords[n].append(coords[n][i] + h * vector_derivative[n])

        i += 1

    return tuple(coords)

"""
Main Script
"""

cid = [0, 1, 3, 6, 0, 4, 5, 8]

def derivative_fnc1(x, y1, y2):
    """
    Method that returns the derivatives of x, y1, and y2 depending on the coordinate.
    Used with an explicit forward euler approximation, this can be used to determine
    what the next values are.

    Returns:
        dx/dx = 1 by definition
        dy1/dx = y2 also by definition
        dy2/dx = -(x + 7)sinx - 5x^2 from simplifying the ode
    """
    y2_x = -(x + 7) * math.sin(x) - (5 * x * y2)
    return 1, y2, y2_x

# Lower bound coordinate = [0, 3rd cid no, 0]
lb_coord = [0, cid[2], 0]
# Lower bound derivative = [1, 5th cid no, 0]
lb_derivative = [1, cid[4], 0]

# Run approximation
x_values, y_values, dy_values = vector_forward_euler(derivative_fnc1, lb_coord, lb_derivative, 15, 0.02)

# Plot values
plot.plot(x_values, y_values)
plot.show()


import numpy, math

# Exercise 1

def backwards_derivative(x_values, y_values):
    """
    Method to return the backwards difference derivative for a set of y values.
    
    y'1 = (y1 - y0) / (x1 - x0)

    The method returns a list one element shorter than the input lists, as the first
    derivative cannot be calculated.
    """
    dy_values = []
    y0 = y_values[0]
    x0 = x_values[0]

    for i in range(1, len(y_values)):
        y1 = y_values[i]
        x1 = x_values[i]
        
        dy_values.append((y1 - y0) / (x1 - x0))

        x0, y0 = x1, y1

    return dy_values

def central_derivative(x_values, y_values):
    """
    Method to return the central difference derivative for a set of y values
    
    y'1 = (y2 - y0) / (x2 - x0)

    The method returns a list two element shorter than the input lists, as the first
    and last derivatives cannot be calculated.
    """
    dy_values = []
    
    y0 = y_values[0]
    x0 = x_values[0]

    for i in range(2, len(y_values)):
        dy_values.append((y_values[i] - y0) / (x_values[i] - x0))

        y0 = y_values[i + 1]
        x0 = x_values[i + 1]

    return dy_values

def create_1_0_array(x, y, z):
    """
    Method to return a 3 dimensional list of any size, with each row either all 1s or
    all 0s. The value of the row is determined by whether the sum of that row's indices 
    are even or odd.
    """
    z_arr = []

    for k in range(z):
        y_arr = []
        for j in range(y):
            y_arr.append([int((j + k) % 2 == 0)] * x)
        z_arr.append(y_arr)

    return z_arr


# Exercise 2

def remove_green(image_array):
    """
    Method to remove all green from an image. This is done by converting the image to a
    two dimensional list, where each element is a tuple containing the values of the RGB
    colours, from 0 to 255.

    The first index is the green value, so this is set to zero.
    A copy of the image is returned so to not affect to input.
    """
    no_green_image_array = image_array.copy()

    for j in range(image_array.shape[0]):
        for i in range(image_array.shape[1]):
            no_green_image_array[j][i][1] = 0
    
    return no_green_image_array

def mirror(image_array):
    """
    Method to return the mirror reflection of the input image.
    A copy of the image is returned so to not affect to input.
    """
    mirror_image_array = image_array.copy()

    for j in range(image_array.shape[0]):
        mirror_image_array[j] = image_array[j][::-1]
    
    return mirror_image_array

def shrink(image_array, n):
    """
    Method to return a smaller version of the input image.
    This is done using a scaling factor, where the method will jump through the image
    array by for each step, rounding to the nearest integer to find the nearest pixel.

    The methods returns a numpy array in the same form as the input.
    """
    small_array = []
    
    for j in range(int(image_array.shape[0] / n)):
        row = []
        for i in range(int(image_array.shape[1] / n)):
            row.append(image_array[int(n * j)][int(n * i)])
        small_array.append(row)

    return numpy.array(small_array)

def vector_derivative(vector, x, y, z, wrt, h):
    """
    Method to return the derivative of a vector with respect to any of the three base 
    vectors using a forward difference approximation.

    If wrt is not 'x', 'y' or 'z' is not calculated (ie it cannot do combinations).
    The method returns a vector in the form of a tuple (i, j, k).
    """
    derivative = [0, 0, 0]

    if h:
        if wrt == 'x':
            p1 = vector(x + h, y, z)
        elif wrt == 'y':
            p1 = vector(x, y + h, z)
        elif wrt == 'z':
            p1 = vector(x, y, z + h)
        else:
            p1 = p0

        p0 = vector(x, y, z)

        for d in (0, 1, 2):
            derivative[d] = (p1[d] - p0[d]) / h

    return tuple(derivative)


def div(vector, x, y, z, h):
    """
    Method to return the divergence of a vector by summing the derivatives of the vector
    in i wrt to x, j wrt to y and k wrt to z.

    Method returns a float, rounded to 5 decimal places.
    """
    i = vector_derivative(vector, x, y, z, 'x', h)[0]
    j = vector_derivative(vector, x, y, z, 'y', h)[1]
    k = vector_derivative(vector, x, y, z, 'z', h)[2]

    return round(i + j + k, 5)

def curl(vector, x, y, z, h):
    """
    Method to return the curl of a vector, returned as a tuple of floats, rounded to
    5 decimal places.
    """
    i = vector_derivative(vector, x, y, z, 'y', h)[2] \
            - vector_derivative(vector, x, y, z, 'z', h)[1]
    j = vector_derivative(vector, x, y, z, 'z', h)[0] \
            - vector_derivative(vector, x, y, z, 'x', h)[2]
    k = vector_derivative(vector, x, y, z, 'x', h)[1] \
            - vector_derivative(vector, x, y, z, 'y', h)[0]

    return (round(i, 5), round(j, 5), round(k, 5))


# Exercise 3

def laplacian(fnc, x, y, z, h):
    """
    Method to calculate the laplacian of a scalar function using a second order central
    difference approximation.
    """
    fx2 = fnc(x + h, y, z)
    fy2 = fnc(x, y + h, z)
    fz2 = fnc(x, y, z + h)

    fx1 = fy1 = fz1 = fnc(x, y, z)

    fx0 = fnc(x - h, y, z)
    fy0 = fnc(x, y - h, z)
    fz0 = fnc(x, y, z - h)

    return ((fx2 - 2 * fx1 + fx0) + \
            (fy2 - 2 * fy1 + fy0) + \
            (fz2 - 2 * fz1 + fz0)) / (2 * h)

def satisfies_laplace(fnc, x_range, y_range, z_range, h):
    """
    Method to return whether a scalar function satisfies the laplace equation. If the
    approximation of the laplace is zero everywhere (rounded to 5 decimal place) then
    True is returned. 
    
    Else the method will stop at the first coordinate where the laplace is not satisfied
    and return False.
    """
    for z in z_range:
        for y in y_range:
            for x in x_range:
                if round(laplacian(fnc, x, y, z, h), 5) != 0:
                    return False
                
    return True

def to_polar(x, y):
    """
    Method to convert a cartesian coordinate to polar
    """
    r = math.hypot(x, y)
    theta = math.atan2(y, x)

    return r, theta


# Exercise 5

def lagrangian_basis_polynomial(i, x, nodes):
    """
    Method to return the polynomial for each y value in the langrangian approximation:
    Lj = product((x - xn)(xj - xn), n = 0...N and n != j) where N = len(nodes)
    """
    polynomial = 1
    for j, (xn, _) in enumerate(nodes):
        if i != j:
            polynomial *= (x - xn) / (nodes[i][0] - xn)
    return polynomial

def lagrangian_approximation(x, nodes):
    """
    Method to return the langrangian approximation at a point x using a set of nodes:
    y = sum((yj * Lj), j = 0...N) where N = len(nodes)
    """
    summation = 0
    for i, (_, y) in enumerate(nodes):
        summation += y * lagrangian_basis_polynomial(i, x, nodes)
    return summation

def ndd(i, nodes):
    """
    Method to return the list of ith forward divided differences for a set of nodes:
    
    (∆^i)y0 = ((∆^i-1)y1 - (∆^i-1)y0) / ((∆^i-1)x1 - (∆^i-1)x0)

    If i = 0, the set of y values is returned. The length of the list returned is equal
    to len(nodes) - i, as the divded differences of the upper nodes can not be calculated.
    """
    if i == 0:
        return list(map(lambda x: x[1], nodes))
    else:
        differences = []
        previous = ndd(i - 1, nodes)

        for j in range(len(previous) - 1):
            differences.append((previous[j + 1] - previous[j]) / (nodes[i + j][0] - nodes[j][0]))

        return differences

def newtonian_polynomial(x, nodes):
    """
    Method to return the newton's polynomial at a point x using a set of nodes:
    y = sum(((∆^n)y * product((x - xi), i = 0...n)), n = 0...N) where N = len(nodes)
    """
    summation = 0

    for i, _ in enumerate(nodes):
        polynomial = 1
        for j in range(i):
            polynomial *= x - nodes[j][0]
        summation += ndd(i, nodes)[0] * polynomial

    return summation

def get_nodes(fnc, x_lb, x_ub, numberof_nodes):
    """
    Method to create a list of nodes to be used for approximating the function fnc.
    The nodes are spaced evenly across the lower and upper x bounds.

    A list of tuples (x, y) is returned.
    """
    coords = []
    x = x_lb
    h = float(x_ub - x_lb) / numberof_nodes

    for _ in range(numberof_nodes):
        coords.append(tuple([x, fnc(x)]))
        x += h

    return coords


# Exercise 6

def adaptive_simpson_approximation(y_values, h):
    """
    Method to return the composite simpson approximation of an integral withequal spacing:
    integral(y, x = a...b) = h/3 * (y0 + yN + 2 * sum(yn, n = 2...N-2, n even) 
        + 4 * sum(yn, n = 1...N-1, n odd)) where N = (b - a) / h
    """
    summation = 0
    for i, y in enumerate(y_values):
        if (i == 0) or (i == len(y_values) - 1):
            summation += y
        elif i % 2 == 0:
            summation += 2 * y
        else:
            summation += 4 * y
    return summation * h / 3


def binomial_coefficient(n, k):
    """
    Method that returns the binomial coefficient n chooses k: (n k) = n! / k!(n - k)!
    Negative and non integer values of n are allowed, k cannot be greater than n if n is
    positive.
    """
    if (n >= 0 and n < k):
        raise ValueError("n must be either negative or greater than k.")
    
    b = float(math.factorial(k)) ** -1

    for i in range(k):
        b *= (n - i)
    
    return b

def binomial_derivative(n, k, nodes, h):
    """
    Method to return an approximation of the nth derivative at a node n using
    the forward equal spacing binomial method.
    
    (yk)^n = sum(((-1)^i * (i k) * y(k + n - i)), i = 0...n)
    """
    summation = 0
    for i in range(n + 1):
        term = nodes[k + (n - i)][1] * binomial_coefficient(n, i)
        summation += ((-1) ** i) * term
    return summation / (h ** n)

def nodal_binomial_derivative(n, nodes):
    """
    Method to return a list of approximations of the nth derivative at each node using
    the forward equal spacing binomial method.

    The list returned will be n elements smaller as the derivatives of the upper nodes
    cannot be calculated.
    """
    derivatives = []
    h = (nodes[-1][0] - nodes[0][0]) / len(nodes)
    for k, (x, _) in enumerate(nodes):
        if k + n < len(nodes) - 1:
            derivatives.append(tuple([x, binomial_derivative(n, k, nodes, h)]))
    
    return derivatives


# Exercise 7

def forward_euler(derivative_fnc, x_lb, y_lb, x_ub, h):
    """
    Method to approximate a function from its derivative using the forward euler method.
    y1 = y0 + hf0
    """
    coords = []
    x, y = x_lb, y_lb
    
    while x < x_ub:
        y = y + (h * derivative_fnc(x, y))
        coords.append(tuple([x, y]))
        x += h
    
    return coords

def runge_kutta_k_values(n, derivative_fnc, coord, h):
    """
    Method to return a list of k values for a runge kutta approximation.
    Method can do any number of runge kutta values, however the most common 
    is 4 (RK4 approximation).
    """
    x, y = coord
    k_values = []
    k = 0
    for i in range(n):
        if i == 0:
            k = h * derivative_fnc(x, y)
        elif i == n - 1:
            k = h * derivative_fnc(x + h, y + k)
        else:
            k = h * derivative_fnc(x + (h / 2), y + (k / 2))
        
        k_values.append(k)

    return k_values

def rk4_approximation(n, derivative_fnc, x_lb, y_lb, x_ub, h):
    """
    Method to return a list of approximate y values using the function derivative and an
    RK4 approximation. y1 = y0 + (h/6)(k1 + 2k2 + 2k3 + k4)
    """
    coords = []
    x, y = x_lb, y_lb
    
    while x < x_ub:
        k = runge_kutta_k_values(4, derivative_fnc, (x, y), h)
        y += (h / 6) * (k[0] + 2 * (k[1] + k[2]) + k[3])
        coords.append(tuple([x, y]))
        x += h
    
    return coords

def vector_forward_euler(derivative_fnc, lb_coord, x_ub, h):
    """
    Method to return a forward euler approximation for a vector. 
    
    As the first element in the coordinates is the domain variable (x, t, etc.) the first
    element in the tuple returned by the function is the derivative of the domain 
    variable with respect to itself and should therefore be equal to 1.
    """
    coords = []
    coord = list(lb_coord)
    
    while coord[0] < x_ub:
        vector_derivative = derivative_fnc(*coord)

        for n in range(len(vector_derivative)):
            coord[n] += h * vector_derivative[n]   
        
        coords.append(tuple(coord))

    return coords


# Exercise 8

def get_compacted_terms(f, g, dx):
    """
    Method to return a set of terms that are easier to work with in numerical calculation.
    Where y'' + f(x)y' + g(x)y = h(x)

    Using finite difference, the terms can be collected into the form:
    (1/dx^2 - f/2dx)y0 + (g - 2/dx^2)y1 + (1/dx^2 + f/2dx)y2 = h
    
    Which are simplified to:
    ay0 + by1 + cy2 = h

    The values a, b and c are returned.
    """
    a = (dx ** -2) - (f / (2 * dx))
    b = g - (2 * (dx ** -2))
    c = (dx ** -2) + (f / (2 * dx))
    
    return a, b, c

def set_boundary_values(type, a, b, c, h, bc, dx, nx):
    """
    Method to return the row values when a boundary condition has to be consiered.

    Using the finite difference equation ay0 + by1 + cy2 = h
    And the robin boundary condition uy' + vy = w.
    These equations can be substituted to remove y-1 from the lower boundary and yn from
    the upper boundary.

    For y-1: (ub/2dx + av)y0 + (ua/2dx + uc/2dx)y1 = uh/2dx + wa
    For yn: (ua/2dx + uc/2dx)y-1 + (ub/2dx - vc)y0 = hu/2dx - wc
    """
    row = numpy.zeros(nx)

    if type == 'lb':
        row[0] = (b * bc[0]) / (2 * dx) + (a * bc[1])
        row[1] = (a + c) * bc[0] / (2 * dx)
        p = (h * bc[0]) / (2 * dx) + (a * bc[2])
    elif type == 'ub':
        row[nx - 2] = (a + c) * bc[0] / (2 * dx)
        row[nx - 1] = (b * bc[0]) / (2 * dx) - (c * bc[1])
        p = (h * bc[0]) / (2 * dx) - (c * bc[2])
    else:
        raise ValueError("Type must be either 'lb' or 'ub'.")
    
    return p, row

def set_inner_values(a, b, c, h, n, nx):
    """
    Method to return the row values when no boundary conditions need to be considered.
    Where ay0 + by1 + cy2 = h.
    """
    row = numpy.zeros(nx)
    
    row[n - 1] = a
    row[n] = b
    row[n + 1] = c

    return h, row

def solve_ode_direct(derivative_fnc, x_lb, x_ub, b_lb, b_ub, nx):
    """
    Method to solve an implicit ode using matrix inversion. The derivative function
    returns three values: f(x), g(x) and h(x) that form the equation:
    
    y'' + f(x)y' + g(x)y = h(x)

    These functions are evaulated at every x value and used to form a matrix of constants
    using the finite difference method:

    (1/dx^2 - f/2dx)y0 + (g - 2/dx^2)y1 + (1/dx^2 + f/2dx)y2 = h

    The boundary conditions b_lb and b_ub are written as robin conditions, where each input
    has three variables u, v and w that form the equation:

    uy' + vy = w --> b_lb = [u1, v1, w1] and b_ub = [u2, v2, w2]

    The method returns a single list of y values corresponding to each x value.
    """
    dx = (x_ub - x_lb) / nx
    
    matrix_of_constants = []
    x_values = []
    p_values = []

    for n in range(nx):
        x = x_lb + (n * dx)
        f, g, h = derivative_fnc(x)
        a, b, c = get_compacted_terms(f, g, dx)

        if n == 0:
            row, p = set_boundary_values('lb', a, b, c, h, b_lb, dx, nx)
        elif n == nx - 1:
            row, p = set_boundary_values('ub', a, b, c, h, b_ub, dx, nx)
        else:
            row, p = set_inner_values(a, b, c, h, n, nx)

        matrix_of_constants.append(row)
        p_values.append(p)
        x_values.append(x)

    y_values = numpy.linalg.inv(matrix_of_constants).dot(p_values)

    return x_values, y_values

def jacobi_approximation(y0, y2, f, g, h, dx):
    """
    Method to return an approximation for y1 using jacobi iteration.

    If ay0 + by1 + cy2 = h, using previous values of y0 and y2, the equation can be
    rearranged to find y1: y1 = (h - cy2 - ay0) / b
    """
    a, b, c =  get_compacted_terms(f, g, dx)
    return (h - (c * y2) - (a * y0)) / b

def jacobi_next_iteration(fnc, x, y, nx, dx):
    """
    Method to run one iteration of the jacobi approximation on previous values to gain a
    better approximation.

    Returns a list of values that apprimxate the function.
    """
    new_y = y.copy()
    
    for n in range(1, nx - 1):
        f, g, h = fnc(x[n])
        new_y[n] = jacobi_approximation(y[n - 1], y[n + 1], f, g, h, dx)
    
    return new_y

def max_difference(x_values, y_values):
    return max(map(lambda x, y: abs(x - y), x_values, y_values))

def solve_ode_jacobi(fnc, x_lb, x_ub, b_lb, b_ub, nx, max_iterations, tolerance):
    """
    Method to solve an ODE iteratively using the jacobi method.
    New values of the function are evaulated by correcting the differences in the
    previous values using the jacobi approximation.

    Returns a list of values when the error between two approximations is smaller than
    the required tolerance, or the approximation reaches a maximum number of iterations.
    """
    x_values = numpy.linspace(x_lb, x_ub, nx, endpoint=True)
    step = (x_ub - x_lb) / nx

    y_values = numpy.zeros(nx)
    previous_y_values = y_values.copy()

    y_values[0] = b_lb
    y_values[nx - 1] = b_ub

    for _ in range(max_iterations):
        y_values = jacobi_next_iteration(fnc, x_values, y_values, nx, step)
        if max_difference(y_values, previous_y_values) < tolerance:
            break
        else:
            previous_y_values = y_values.copy()

    return x_values, y_values


# Exercise 9

def is_stable(dt, dx, a):
    """
    Method that returns if the time and space steps of a poisson approximation is valid.
    """
    return a * dt * (dx ** -2) < 0.5

def central_point_temperature(tx0y1, tx1y0, tx1y1, tx1y2, tx2y1, ds, dt, a):
    """
    Method to approximate the new temperature at a point by a heat flux due to
    temperature differences in both the x and y directions, using the finite difference
    method.
    """
    heat_flux_x = a * (tx2y1 - (2 * tx1y1) + tx0y1) / (ds ** 2)
    heat_flux_y = a * (tx1y2 - (2 * tx1y1) + tx1y0) / (ds ** 2)

    return tx1y1 + dt * (heat_flux_x + heat_flux_y)

def get_temperature_frame(phi, fnc, p, ds, dt, nx):
    """
    Method that returns the temperature distribution at a single moment in time.
    """
    frame = []
    
    for j in range(nx):
        row = []
        for i in range(nx):
            row.append(fnc(phi[p - 1], i, j, ds, dt, nx))
        
        frame.append(row)

    return frame
        

def temperature_distribution(phi_fnc, initial_temp, x_lb, y_lb, x_ub, y_ub, t_end, nt, nx):
    """
    Method that returns the temperature distribution varying with time in the form of a
    three dimensional list.
    """
    ds = (x_ub - x_lb) / nx
    dt = t_end / nt

    phi = [initial_temp]
    
    for p in range(1, nt):
        phi.append(get_temperature_frame(phi, phi_fnc, p, ds, dt, nx))

    return phi

def example_phi_function(phi, i, j, ds, dt, nx):
    """
    An example method that defines the temperature distribution of an object.
    All boundary definitions are given at an edge and are described by a function that
    returns the temperature based on location or the type of condition set.

    Away from the boundaries conditions, a finite difference approximation based on the
    poisson equation is used.
    """
    a = 1 # Thermal diffusivity
    if i == 0: # Example constant temperature edge
        return 0
    elif j == nx - 1: # Example adiabatic edge
        return adiabatic_boundary(phi, i, j, ds, dt, a)
    elif j == 0: # Example varying temperature edge
        return 273.15 * math.sin(math.pi * i / nx)
    elif i == nx - 1: # Example robin bc edge
        return mixed_boundary(phi, i, j, ds, dt, a)
    else:
        return central_point_temperature(
            phi[i-1][j], phi[i][j-1], phi[i][j], phi[i][j+1], phi[i+1][j], ds, dt, a)

def adiabatic_boundary(phi, i, j, ds, dt, a):
    """
    Method that returns the temperature at a point along an edge with an adiabatic
    boundary condition: y' = 0 --> (yi1j2 - yi1j0) / 2dx = 0 --> yi1j2 = yij0

    This particular case is for a vertical edge, therefore this does not apply to the
    x direction. The rest of the calculation is carried out using the poisson equation.
    """
    phi_i1_j2 = phi[i][j-1]
    
    return central_point_temperature(
        phi[i-1][j], phi[i][j-1], phi[i][j], phi_i1_j2, phi[i+1][j], ds, dt, a)

def mixed_boundary(phi, i, j, ds, dt, a):
    """
    Method that returns the temperature at a point along an edge with a mixed
    boundary condition:
    uy' + vy = w --> u * (yi2j1 - yi2j1) / 2 dx = w - v * yi1j1 
                 --> yi2j1 = yi0j1 + (2 dx / u) * (w - v * yi1j1)

    Where u, v and w are constants of the robin boundary condition.

    This particular case is for a horizontal edge, therefore this does not apply to the
    y direction. The rest of the calculation is carried out using the poisson equation.
    """
    u, v, w = 1, 1, 1
    phi_i2_j1 = phi[i-1][j] + (2 * ds / u) * (w - v * phi[i][j])
    
    return central_point_temperature(
        phi[i-1][j], phi[i][j-1], phi[i][j], phi[i][j+1], phi_i2_j1, ds, dt, a)


# Coursework

def get_cfl(c, dt, dx):
    """
    Calculates the CFL of a finite difference method.
    CFL = c * (dt / dx)
    """
    return c * dt / dx

def central_wave_amplitude(u, i, p, cfl):
    """
    Explicit finite difference method for calculating the amplitude
    of a point i along a string (1D) given that p > 1:

    u(i,p) = (2u(i,p-1) - u(i,p-2))
        + cfl^2 (u(i-1,p-1) - 2u(i,p-1) - u(i+1,p-1))
    """
    finite_diff_time = (2 * u[p - 1][i]) - u[p - 2][i]
    finite_diff_x = u[p - 1][i - 1] - (2 * u[p - 1][i]) + u[p - 1][i + 1]
    
    return finite_diff_time + (cfl ** 2) * finite_diff_x

def initial_wave_amplitude(u, i, dt, v0):
    """
    Explicit method for calculating the amplitude of a point i
    along a string (1D) at p = 1:

    u(i,p) = u(i,p-1) + dt * u'(i,p-1)
    """
    return u[0][i] + dt * v0[i]

def wave_frame(u, v0, p, bc, cfl, nx, dt):
    """
    Method to return the string amplitudes at a given frame.
    
    If a gridpoint is in the list of boundary conditions, the value of that
    boundary is added to the list. Therefore all boundary conditions must be
    Dirchlet type.
    
    If p = 1, use initial_wave_amplitude(), else use normal explicit finite
    difference method central_wave_amplitude().
    """
    frame = []
    
    for i in range(nx):
        if i in bc:
            frame.append(0) # all boundary values are zero
        elif p == 1:
            frame.append(initial_wave_amplitude(u, i, dt, v0))
        else:
            frame.append(central_wave_amplitude(u, i, p, cfl))

    return frame

def get_wave_oscillation(length, t_end, u0, v0, bc, c, nx, nt):
    """
    Method to return a list of frames of the string oscillating.
    As each frame is a set of y values, the method returns a list of lists.
    """
    u = [u0]
    dt = t_end / nt
    dx = length / nx

    cfl = get_cfl(c, dt, dx)

    for p in range(1, nt):
        u.append(wave_frame(u, v0, p, bc, cfl, nx, dt))
    
    return u
from matplotlib import pyplot as plot
import math as m

# Task A

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

# Task B

def binomial_coefficient(n, k):
    return m.factorial(n) / (m.factorial(k) * m.factorial(n - k))

def node_binomial_derivative(n, k, y_values, h):
    summation = 0
    for i in range(n + 1):
        term = y_values[k + (n - i)] * binomial_coefficient(n, i)
        summation += ((-1) ** i) * term
    return summation / (h ** n)

def binomial_derivative(n, nodes):
    derivatives = []
    h = (nodes[-1][0] - nodes[0][0]) / len(nodes)
    y_values = list(map(lambda x: x[1], nodes))
    for k, (x, _) in enumerate(nodes):
        if k + n < len(nodes) - 1:
            derivatives.append(tuple([x, node_binomial_derivative(n, k, y_values, h)]))
    
    return derivatives

# Task C

def get_nodes(fnc, x_lb, x_ub, numberof_nodes):
    coords = []
    x = x_lb
    h = float(x_ub - x_lb) / numberof_nodes

    for _ in range(numberof_nodes):
        coords.append(tuple([x, fnc(x)]))
        x += h

    return coords

# Approximation using 100 nodes
sinx_nodes = get_nodes(m.sin, 0, m.pi, 100)
sinx_5th_derivative = binomial_derivative(5, sinx_nodes)

plot.plot(
    list(map(lambda x: x[0], sinx_nodes)), list(map(lambda x: x[1], sinx_nodes)),
    list(map(lambda x: x[0], sinx_5th_derivative)), list(map(lambda x: x[1], sinx_5th_derivative))
)

plot.show()
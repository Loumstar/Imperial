import math as m
from matplotlib import pyplot as plot

def lagrangian_basis_polynomial(i, x, nodes):
    polynomial = 1
    for j, (xn, _) in enumerate(nodes):
        if i != j:
            polynomial *= (x - xn) / (nodes[i][0] - xn)
    return polynomial

def lagrangian_approximation(x, nodes):
    summation = 0
    for i, (_, y) in enumerate(nodes):
        summation += y * lagrangian_basis_polynomial(i, x, nodes)
    return summation

def ndd(i, nodes):
    if i == 0:
        return list(map(lambda x: x[1], nodes))
    else:
        differences = []
        previous = ndd(i - 1, nodes)
        
        for j in range(len(previous) - 1):
            differences.append((previous[j + 1] - previous[j]) / (nodes[i + j][0] - nodes[j][0]))
    
        return differences

def newtonian_polynomial(x, nodes):
    summation = 0
    
    for i, _ in enumerate(nodes):
        polynomial = 1
        for j in range(i):
            polynomial *= x - nodes[j][0]
        summation += ndd(i, nodes)[0] * polynomial

    return summation

def get_nodes(fnc, x_lb, x_ub, numberof_nodes):
    coords = []
    x = x_lb
    h = float(x_ub - x_lb) / numberof_nodes

    for _ in range(numberof_nodes):
        coords.append(tuple([x, fnc(x)]))
        x += h

    return coords

def get_approximation_distribution(nodes, x_lb, x_ub, h):
    coords = []
    x = x_lb
    while x < x_ub:
        coords.append(tuple([x, newtonian_polynomial(x, nodes)]))
        x += h
    return coords

sine_nodes = get_nodes(m.sin, 3, 5, 4)
sine_coords = get_nodes(m.sin, 0, 8, 100)

sine_approximation = get_approximation_distribution(sine_nodes, 0, 8, 0.01)

plot.plot(
    list(map(lambda x: x[0], sine_approximation)), list(map(lambda x: x[1], sine_approximation)),
    list(map(lambda x: x[0], sine_coords)), list(map(lambda x: x[1], sine_coords))
)

plot.show()
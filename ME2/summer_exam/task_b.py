import numpy
from matplotlib import pyplot as plot

cid_values = [0, 1, 3, 6, 0, 4, 5, 8]
cid_x_values = numpy.arange(0.8, step=0.1)

# Create a large array of numbers to be used to plot the approximation
lagrange_x_values = numpy.linspace(0, 0.7, num=100, endpoint=True)

lagrange_y_values = []

for x in lagrange_x_values:
    
    summation = 0
    
    for i, y in enumerate(cid_values):
        
        # Calculate the lagrangian basis polynomial for each node
        # where Lj = product((x - xn)(xj - xn), n = 0...N and n != j)
        polynomial = 1
        for j, xn in enumerate(cid_x_values):
            if i != j:
                polynomial *= (x - xn) / (cid_x_values[i] - xn)
        
        # Multiply each basis polynomial with the y value of the node and
        # sum these to gain the approximation:
        # y = sum((yj * Lj), j = 0...N) where N = len(nodes)
        summation += y * polynomial

    # Add to list of approximations
    lagrange_y_values.append(summation)

# Plot cid numbers
plot.scatter(cid_x_values, cid_values)
# Plot lagrange approximation
plot.plot(lagrange_x_values, lagrange_y_values, 'r')

plot.show()
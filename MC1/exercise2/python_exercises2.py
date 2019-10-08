import numpy as np, matplotlib.pyplot as plot

print("--- Task A")
#1.
A = [x for x in range(11, 21)]
B = [y for y in range(21, 31)]
#2.
A[4] = A[2] + A[3]
B[5] *= 2
#3.
A[0], A[len(A)-1] = (A[len(A)-1], A[0])
#4.
i, j = (3, 5)
B[j], A[i] = (A[i], B[j])
print(A[0])
print("--- Task B")
#1.
N = 100
A = [x for x in range(1, N+1)]
#2.
B = [y**2 for y in A]
#3.
C = [A[i]+B[i] for i in range(len(A))]
print(C)
print(C[25])
print("--- Task C")
#1.
A = list(
    map(lambda x, i: x if (i+1) % 3 != 0 else 0, 
    A, 
    [i for i in range(len(A))])
)
#2.
D = [A[len(A)-j-1] for j in range(len(A))]
print(D[9])
print("--- Task D")
x_plot = [(x+1)*5 for x in range(20)]
y_plot = [np.log10(x) for x in x_plot]

plot.plot(x_plot, y_plot)
plot.ylabel("Log10(x)")
plot.xlabel("x")

plot.show()
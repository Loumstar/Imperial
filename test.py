zero_martrix = [
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0]
]

def zero_matrix(n, m):
    mat = []
    for _ in range(n):
        row = [0] * m
        mat += [row]
    return mat

def multiply_matrices(a, b):
    #find number of sums that have to be made to find each element of c
    if len(a[0]) == len(b):
        d = len(b)
    else:
        return None
    
    c = []
    m, n = len(a), len(b[0])

    for i in range(m):
        row = []
        for j in range(n):
            s = 0
            for term in range(d):
                s += a[i][term] * b[term][j]
            row += [s]
        c += [row]
    
    return c

a1 = [
    [1, 2],
    [3, 4],
    [5, 6]
]

b1 = [
    [9, 8, 7],
    [6, 5, 4]
]

print(multiply_matrices(a1, b1))
            
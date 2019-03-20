import re

print("--- Task A")

def fib(a, b): 
    return a + b

def nth_fibonacci(n):
    a = b = 1
    for _ in range(n - 1):
        a, b = b, fib(a, b)
    return a

def recursive_fibonacci(n):
    if n < 2:
        return 1
    else:
        return recursive_fibonacci(n-1) + recursive_fibonacci(n-2)

print(recursive_fibonacci(4))
#print(nth_fibonacci(22))

print("--- Task B")

m = [
    [24, 67, 81],
    [10, 3, 28],
    [63, 55, 17]
]

def swap_elements(a, b, m):
    #where a and b are tuples of ij index notation
    m[a[0]][a[1]], m[b[0]][b[1]] = m[b[0]][b[1]], m[a[0]][a[1]]
    return m

val67 = (0, 1)
val17 = (2, 2)
val81 = (0, 2)
val63 = (2, 0)

swap_elements(val67, val17, m)
swap_elements(val81, val63, m)

m[1][1] = m[0][0] + m[1][2]

print(m[2][2])

print("--- Task C")

def create_flag(n):
    m = []
    for j in range(n):
        sub_m = [0] * n
        sub_m[j] = 1
        m.append(sub_m)
    return m

print("--- Task D")

with open("./exercise7/MatA.txt", 'r') as f:
    mat_a_lst = [int(x) for x in re.findall('\d+', f.read())]

with open("./exercise7/MatB.txt", 'r') as f:
    mat_b_lst = [int(x) for x in re.findall('\d+', f.read())]

def read_matrix(lst, n, m):
    matrix = []
    for j in range(m):
        i = j * n
        matrix.append(lst[i:i+n])
    return matrix


mat_a = read_matrix(mat_a_lst, 60, 60)
mat_b = read_matrix(mat_b_lst, 60, 60)

print(mat_b[25][36])

def sum_matrix(a, b):
    m = []
    for j in range(len(a)):
        m.append(list(map(lambda x, y: x + y, a[j], b[j])))
    return m

print(sum_matrix(mat_a, mat_b)[24][19])

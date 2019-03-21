def lattice_paths(n, m, *position):
    """
    Recursive function to determine how many unique paths can be
    taken to reach the bottom-right corner, from the top-left corner,
    of a m x n square
    """
    position = (0, 0) if position == () else position[0]
    
    total_d = 0
    total_r = 0

    if position == (n, m):
        print(position)
        return 1

    #Right option
    if position[0] != n:
        position_r = position[0] + 1, position[1]
        total_r = lattice_paths(n, m, position_r)

    #Down option
    if position[1] != m:
        position_d = position[0], position[1] + 1
        total_d = lattice_paths(n, m, position_d)
        
    return total_d + total_r

def new_lattice_paths(n, m):
    """
    Second version of the script, which still uses recursion
    but skips a few steps to increase its speed
    """
    largest_dim = n if n > m else m
    smallest_dim = m if m > n else n
    
    s = sum_paths(largest_dim, smallest_dim)
    s += largest_dim - smallest_dim
    
    return s

def sum_paths(n, r):
    print(r)
    s = 0
    if r <= 1:
        return n + 1
    else:
        for j in range(n+1):
            s += sum_paths(j, r-1)
    return s

def binomial_lattice_paths(n, m):
    """
    Third and final version of the script, which takes
    a whole new approach, using binomial theorem.

    Part of the lattice is ignored at first so that it is a perfect square of size n x n or m x m, whichever is smaller.
    The number of paths possible given the first down movement is at the k-th column from the right is equal to
    the sum of paths possible given the first down movement is at (k-1)-th column from right
    for a n-1 x n-1 square, plus a n-2 x n-2 square and so on til square has no size.

                <--- k
    +----+----+----+                 +----+----+----+                 +----+----+               +----+----+
    |    |    |    |                                |                           |                         |
    |----|----|----|       k = 0                    |       k = 1               +----+                    |
    |    |    |    | n    ------>                   | n    ------>                   |      +             +----+     +   ...
    |----|----|----|                                |                                | x                       |
    |    |    |    |                                |                                |                         | x
    +----+----+----+                                |                                |                         |
            n                                                               x = n-1                 x = n-2
    
               +----+             +----+
     k = 1          |                  |
    ------->        | n-1    +         | n-2    +   ... which is equivalent to the sum of smaller squares with k = 0
                    |                  |
                    |


    In the example above that for k = 1, for each value of x, the sum of the paths can be broken down into
    the sum of paths of smaller squares with k = k - 1

    ie:
    sum_paths(n, k) = sum_paths(n-1, k-1) + sum_paths(n-2, k-1) + ... + sum_paths(0, k-1)
    
    This is equivalent to summing the binomial(x, k-1), where x = [0...n-1].

    However, this value already exists in pascals triangle as binomial(n+k, k+1), 
    which is what this script finds for each value of k. These values are summed and then returned.

    The number of extra paths due to the lattice being retangular is equal to the difference in the dimensions.
    """    
    largest_dim = n if n > m else m
    smallest_dim = m if m > n else n
    s, k = 0, 0
    while k <= smallest_dim:
        s += binomial(smallest_dim + k - 1, k)
        k += 1
    s += largest_dim - smallest_dim
    return s

def fact(a):
    p = 1
    for i in range(a):
        p *= i + 1
    return p

def binomial(n, k):
    return int(fact(n) / (fact(k) * fact(n-k)))

print("binomial", binomial_lattice_paths(10, 10))
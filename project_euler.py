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


#print(lattice_paths(10, 10))

def new_lattice_paths(n, m):
    largest_col = n if n > m else m
    smallest_col = m if m > n else n
    
    s = sum_paths(largest_col, smallest_col)
    s += largest_col - smallest_col
    
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

print(new_lattice_paths(3, 3))
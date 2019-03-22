from random import random as r
import math as m
from matplotlib import pyplot as plot

# Question 1
def matrix_product(a, b):
    """
    Method to determine the multiplcation of two matrices a and b.
    If the matrices are not compantible, None is returned.
    If they are compatible, each element in row a is multiplied by 
    the element at that index in column b. These are summed to give 
    the value of the element of the end-result matrix c, at that 
    row and column location.
    """
    #checks that the number of columns in a is equal to the number of rows in b
    if len(a[0]) != len(b):
        #if they are not equal, return None and exit
        return None
    #determine the dimensions of the new matrix
    #n number of rows, m number of columns
    n, m = len(a), len(b[0])
    #the resultant matrix is created as an empty list
    c = []
    #for each row of c
    for j in range(n):
        row = []
        #for each element in that row. Alternatively, for each column in b
        for i in range(m):
            #multiply each element in row a with each element in column b that is at the same index.
            #find the sum, and append that sum to the row list.
            row.append(sum(list(
                map(lambda x, y: x*y[i], a[j], b)
            )))
        #append each row to c, to give a list of lists.
        c.append(row)
    return c

def transpose(m):
    """
    Method to return the transpose of a matrix.
    The matrix must be square. The columns are 
    found and then appended as a row to a new matrix r.
    """
    #create an empty list for resultant matrix
    r = []
    #for each column in m
    for j in range(len(m[0])):
        #append, to r, a list of all elements that are at index j of each row
        #i.e. append a list of the elements in one column to r as a row.
        r.append(list(map(lambda x: x[j], m)))
    return r

def multiply_matrix_by_transpose(m):
    #find the transpose
    t = transpose(m)
    #find the product of the original matrix, by its transpose
    return matrix_product(m, t)

a = [
    [1, 2, 3, 4],
    [4, 5, 6, 7],
    [8, 9, 8, 7],
]

b = [
    [1],
    [2],
    [3],
    [4]
]

print(matrix_product(a, b))

# Question 2
def create_brick(lb, ub, length):
    """
    Method to create a data point random between bounds that form a square.
    These data points are then companred to length to ensure no part of the
    brick falls outside the square boundary, else new data points will be
    generated.
    """
    #lb and ub are tuples of coordinates, indicating the lower 
    #and upper bounds for where the brick can exist
    #method picks a random point in the range of the upper and lower bounds
    #and adds that to the lower bound for each coordinate
    x = (r() * (ub[0] - lb[0])) + ub[0]
    y = (r() * (ub[1] - lb[1])) + ub[1]

    #if the length causes the brick to fall outside the boundary 
    #pick two new random points
    while x - (length / 2) < lb[0] or x + (length / 2) > ub[0]:
        x = (r() * (ub[0] - lb[0])) + lb[0]
        y = (r() * (ub[1] - lb[1])) + lb[1]        

    return x, y


def draw_brick(brick_coords, length):
    """
    Method to create a set of data points that will create the image of
    the brick in pyplot.
    """
    #line data points starts as a tuple of empty lists for x and y coords
    line = ([], [])
    for i in range(11):
        #add ten data points equidistant along the length of the brick
        #note brick coord located at the centre of the brick
        line[0].append(brick_coords[0] + (i * length / 10) - (length / 2))
        #add same y coord
        line[1].append(brick_coords[1])
    return line

def draw_shooting_line(theta, y_coord):
    """
    Method to create a set of data points that will create the image of
    the shooting line in pyplot.
    """
    #line data points starts as a tuple of empty lists for x and y coords
    line = ([], [])
    for i in range(10):
        #plot y coordinates that are between the origin and the brick
        line[1].append(y_coord / (i + 1))
        #plot x coordinates with respect to y values: y = x tan(a)
        line[0].append(y_coord / ((i + 1) * m.tan(theta)))
    return line

def deg(theta):
    #function to convert radians to degrees
    return theta * 180 / m.pi

def find_shooting_angles(brick_coords, length):
    """
    Method to find the upper and lower bounds of the shooting
    angle that will allow the brick to be hit.
    """
    #theta is larger when x is smaller, therefore max theta is
    #the angle from the origin to the closest end of the brick, 
    #with respect to its height

    #min theta is the opposite, the angle from origin to the furthest
    #end of the brick, with resect to its height.
    min_theta = m.atan(brick_coords[1] / (brick_coords[0] + (length / 2)))
    max_theta = m.atan(brick_coords[1] / (brick_coords[0] - (length / 2)))

    return min_theta, max_theta

def Brick(lb, ub, length):
    #create brick coordinates
    brick = create_brick(lb, ub, length)

    #determine shooting angles and the midpoint angle
    shooting_angles = find_shooting_angles(brick, length)
    mid_angle = sum(shooting_angles) / 2

    #draw the brick and the line of action
    brick_line = draw_brick(brick, length)
    shooting_line = draw_shooting_line(mid_angle, brick[1])

    #plot the lines, the zeroth index refers to the x coordinate and
    #the first index refers to the y coordinate
    plot.plot(
        brick_line[0], brick_line[1], 'r',
        shooting_line[0], shooting_line[1], 'b',
        brick[0], brick[1], 'gx'
    )
    plot.show()

Brick((0, 0), (2, 2), 0.2)

# Question 3
with open("./exam/FTSE100.txt", "r") as f:
    #create a list for each line
    ftse100_list = f.read().split("\n")

def companies(ftse):
    """
    Method to sort the list of companies not yet organised into tuples.
    """
    #a list that will contain tuples of information for each company
    companies = []
    #for every third line in the list, as each company takes up 
    #three lines of information
    for c in range(int(len(ftse) / 3)):
        #add the three lines as separate elements in a tuple 
        #to the compant list
        companies.append(
            (ftse[(3 * c)], ftse[(3 * c) + 1], ftse[(3 * c) + 2])
        )
    return companies

def sort_by_share_price(c):
    #assume first company has the largest share price
    c_sorted = [c[0]]
    #for all companies in list
    for x in c[1:]:
        #assume they are the smallest
        smallest = True
        #for each company in the sorted list
        for i,y in enumerate(c_sorted):
            #if the company being tested has a higher share price
            if float(x[1]) > float(y[1]):
                smallest = False
                #the test company will be moved to that index
                #shifting the rest of the list one forward
                c_sorted.insert(i, x)
                break
        if smallest:
            #simply add the company on the end of the list
            c_sorted.append(x)
    
    return c_sorted

def total_ftse_value(companies):
    #initial sum is zero
    s = 0.0
    for c in companies:
        #multiply share price by number of share available
        s += float(c[1]) * float(c[2])
    #divide by one billion
    return s * (10 ** -9)

def LSE(ftse):
    #organise the company information
    companies_list = companies(ftse)
    #sort by share price
    companies_sorted = sort_by_share_price(companies_list)
    #determine the value of the FTSE 100
    ftse_value = total_ftse_value(companies_list)
    
    #return sorted list and ftse value as a tuple
    return companies_sorted, ftse_value

print(LSE(ftse100_list))
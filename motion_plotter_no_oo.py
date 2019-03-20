from random import uniform as rand
from matplotlib import pyplot as plt
import math

def draw_boundary(radius):
    #create a tuple of lists, which will take the data points of the boundary
    boundary = ([], [])

    for a in range(361):
        #add a data point for each degree of the cirle that makes up the boundary
        boundary[0].append(radius * math.sin(a * math.pi / 180))
        boundary[1].append(radius * math.cos(a * math.pi / 180))

    return boundary
        
def random_data_point(radius):
    #creates random data points in a square betweem (-radius, -radius) and (radius, radius) 
    return (rand(-1, 1) * radius), (rand(-1, 1) * radius)

def create_ant(radius):
    initial_x, initial_y = random_data_point(radius)

    #checks the values are within the boundary, if not create some new ones
    while (initial_x ** 2) + (initial_y ** 2) > (radius ** 2):
        initial_x, initial_y = random_data_point(radius)
    
    #ant is a tuple of lists containing the data points, starting with its initial position.
    ant = ([initial_x], [initial_y])

    return ant

def move(jump, radius, ant):
    dx = rand(-1, 1) * jump
    dy = rand(-1, 1) * jump

    #find the latest ant position, as the last item in the lists of the ant tuple
    ant_position_x = ant[0][len(ant[0])-1]
    ant_position_y = ant[1][len(ant[1])-1]
    
    new_x = ant_position_x + dx
    new_y = ant_position_y + dy
    
    #if the ant jumps only if it doesn't all out side the boundary
    if (new_x ** 2) + (new_y ** 2) <= (radius ** 2):
        ant[0].append(new_x)
        ant[1].append(new_y)

    return ant

#boundary has radius 1
boundary = draw_boundary(1)
#ant created so its initial position is within radius 1
ant = create_ant(1)

for _ in range(150):
    ant = move(0.2, 1, ant)

plt.plot(
    boundary[0], boundary[1], 'b',
    ant[0], ant[1], 'r',
    ant[0][0], ant[1][0], 'gx'
)

plt.show()

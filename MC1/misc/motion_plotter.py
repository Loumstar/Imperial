from random import uniform as r
from matplotlib import pyplot as plt
import math

class Ant:
    def __init__(self, jump, radius):
        self.jump = jump
        self.radius = radius
        self.boundary = {'x': [], 'y': []}
        self.ant = {'x': [], 'y': []}
        self.ant_position = (None, None)
        self.initial_position = (None, None)

    def draw_boundary(self):
        for a in range(361):
            self.boundary['x'].append(self.radius * math.sin(a * math.pi / 180))
            self.boundary['y'].append(self.radius * math.cos(a * math.pi / 180))

        return self.boundary
        
    def random_data_point(self):
        return (r(-1, 1) * self.radius), (r(-1, 1) * self.radius)

    def create_ant(self):
        initial_x, initial_y = self.random_data_point()
        
        while (initial_x ** 2) + (initial_y ** 2) > (self.radius ** 2):
            initial_x, initial_y = self.random_data_point()

        self.ant_position = (initial_x, initial_y)
        self.initial_position = (initial_x, initial_y)
        
        self.ant['x'].append(initial_x)
        self.ant['y'].append(initial_y)

        return self.ant_position

    def move(self):
        dx = r(-1, 1) * self.jump
        dy = r(-1, 1) * self.jump
        
        new_x = self.ant_position[0] + dx
        new_y = self.ant_position[1] + dy
        
        if not (new_x ** 2) + (new_y ** 2) > (self.radius ** 2):
            self.ant['x'].append(new_x)
            self.ant['y'].append(new_y)

            self.ant_position = (new_x, new_y)

        return self.ant_position

    def simulate(self, n):
        if self.ant_position != (None, None) and self.boundary != {'x': [], 'y': []}:
            for _ in range(n):
                self.move()
        else:
            print("Ant position or boundary have not been initialised")

    def plot(self):
        plt.plot(
            self.boundary['x'], self.boundary['y'], 'b',
            self.ant['x'], self.ant['y'], 'r',
            [self.initial_position[0]], [self.initial_position[1]], 'gx'
        )
        plt.show()


test_ant = Ant(0.2, 1)

test_ant.draw_boundary()
test_ant.create_ant()
test_ant.simulate(1000)
test_ant.plot()
from matplotlib import pyplot as plot
from random import randint as r
import math
import re

print("--- Task A")

def factorial(x):
    f = 1
    for i in range(x):
        f *= i + 1
    return f

print(factorial(5))

print("--- Task B")

def series_expansion(r, n):
    data_points = []
    for x in r:
        y = 0
        for p in range(n+1):
            y += (x ** p) / factorial(p)
        
        data_points.append(y)
    
    return data_points

example_range = list(map(lambda x: x / 100, range(500)))

plot.plot(
    example_range, series_expansion(example_range, 2),
    example_range, series_expansion(example_range, 6),
    example_range, series_expansion(example_range, 10),
    example_range, series_expansion(example_range, 14)
)

print(series_expansion(example_range, 18)[example_range.index(4.8)])

print("--- Task C")

def sort_ascending(l):
    sorted_list = [l[0]]
    for x in l[1:]:
        smallest = True
        for i,y in enumerate(sorted_list):
            if x > y:
                smallest = False
                sorted_list.insert(i, x)
                break
        if smallest:
            sorted_list.append(x)
    
    return sorted_list

with open('./Set.txt') as f:
    blackboard_numbers = [int(x) for x in re.findall('\d+', f.read())]

sorted_bb_numbers = sort_ascending(blackboard_numbers)
print(list(reversed(sorted_bb_numbers))[24])

def test_sort_ascending(l):
    v = l[0]
    for x in l:
        if x > v:
            return False
        v = x
    return True

random_numbers = [r(1, 100) for _ in range(50)]

print(
    test_sort_ascending(
        sort_ascending(random_numbers)
    )
)

print("--- Task D")

def free_fall(x, a, v, g):
    return (x * math.tan(a)) - (((x / (v * math.cos(a))) ** 2) * g / 2)

def draw_target(target):
    data_points = {'x': [], 'y': []}
    
    for a in range(361):
        data_points['x'].append(target['r'] * (1 + math.sin(a * math.pi / 180)))
        data_points['y'].append(target['r'] * (1 + math.cos(a * math.pi / 180)))
    
    return data_points

def shoot_bullet(theta, v, target, step):
    data_points = {'x': [], 'y': []}
    target_data_points = draw_target(target)
    
    x = 0
    a = theta * math.pi / 180
    bullseye = False

    while x < target['x'] + target['r']:
        y = free_fall(x, a, v, 9.81)
        
        data_points['x'].append(x)
        data_points['y'].append(y)

        if ((x - target['x']) ** 2) + ((y - target['y']) ** 2) < target['r'] ** 2:
            bullseye = True
            break
    
        x += step

    plot.clf()

    plot.plot(
        data_points['x'], data_points['y'],
        target_data_points['x'], target_data_points['y']
    )

    plot.show()

    return bullseye

example_target = {
    'r': 0.01,
    'x': 0.3,
    'y': 0.7
}

print(shoot_bullet(63, 21, example_target, 0.01))
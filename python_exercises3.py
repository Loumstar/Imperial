from random import randint, random
from matplotlib import pyplot as plot
from math import pi

print("--- Task A")
def loop_phrase(n):
    for x in range(n):
        print("I need to comment my scripts. Comments are marked too!")
    
def sum_x2(n):
    sum = 0
    for x in range(n):
        sum += (x+1) ** 2
    return sum

def sum_dice(n):
    sum = 0
    for x in range(n):
        sum += randint(1, 6)
    return sum

def factorial(n):
    x = n
    product = 1
    while x > 0:
        product *= x
        x -= 1
    return product

print("Sum x^2", sum_x2(3))
print("Sum dice", sum_dice(4))
print("Factorial", factorial(13))

print("--- Task B")

class PlayDice:
    def __init__(self):
        self.player_a = input("Type Player A\'s name: ")
        self.player_b = input("Type Player B\'s name: ")
        self.total_a = 0
        self.total_b = 0

    def play(self):
        a = randint(1, 6)
        b = randint(1, 6)
        
        print("%s: %i, %s: %i" % (self.player_a, a, self.player_b, b))
        if(a > b):
            print("%s wins %s" % (self.player_a, self.player_b))
            self.total_a += 1
        else:
            print("%s wins %s" % (self.player_b, self.player_a))
            self.total_b += 1

    def report(self):
        print("%s won %i games, while %s won %i games." % (self.player_a, self.total_a, self.player_b, self.total_b))

def rock_paper_scissors():
    rps = ["r", "p", "s"]
    user = rps.index(
        input("Enter Rock (R), Paper (P) or Scissors (S): ").lower()
    )
    computer = randint(1, 3)
    if computer == user:
        print("Draw!")
    elif computer > user or (computer == 1 and user == 3):
        print("Computer wins!")
    else:
        print("User wins!")

game = PlayDice()

game.play()
game.play()
game.report()

rock_paper_scissors()

print("--- Task C")

def calculate_pi(n):
    """
    The method is to choose random numbers between 0 and 1 in both the x and y plane, representing a square.
    
    Each x and y coord is saved to the square list unless the line from origin to the point is smaller 
    than the diameter of a circle diameter 1, it is saved to a separate list.
    if(x**2 + y**2 > r**2) append to circle list
    
    The random numbers land in or out of the circle depending on the ratio of the areas of the circle and square.
    Therefore the ratio of data points in the circle and the square will approximate pi/4.
    """
    points_in_circle = 0
    x_points = {"square":[], "circle":[]}
    y_points = {"square":[], "circle":[]}
    i = 0
    while i < 10**n:
        x = random()-0.5
        y = random()-0.5

        if x**2 + y**2 < 0.25:
            points_in_circle += 1
            x_points["circle"].append(x)
            y_points["circle"].append(y)
        else:
            x_points["square"].append(x)
            y_points["square"].append(y)

        i += 1

    pi_approximation = 4 * points_in_circle / (10**n)

    print("Pi approximation: %f, Error: %f" % (pi_approximation, pi_approximation - pi))

    plot.plot(x_points["circle"], y_points["circle"], 'rx', x_points["square"], y_points["square"], 'bx')
    plot.ylabel("y")
    plot.xlabel("x")

    plot.show()

calculate_pi(4)

print("--- Task D")

def list_primes_to(n):
    primes = [2]
    i = 3
    while i < n:
        is_prime = True
        for x in primes:
            if i % x == 0:
                is_prime = False
                break
        if is_prime:
            primes.append(i)
        i += 1
    return "Found %i primes!" % len(primes)

print(list_primes_to(2680))

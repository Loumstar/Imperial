import time as t

def force(G, M, m, r):
    return (G * M * m) / (r ** 2)

G = -6.67 * (10 ** -11)

mass_1 = 1.989 * (10 ** 30)
mass_2 = 1
displacement = 150 * (10 ** 9)
velocity = 0
time = 0

delta_t = 100

start = t.time()

while displacement > 0:
    print(round(displacement, 2), round(velocity, 2), round(time, 2))

    acceleration = force(G, mass_1, mass_2, displacement) * ((1/mass_1) + (1/mass_2))
    velocity += acceleration * delta_t
    displacement += velocity * delta_t
    
    time += delta_t

days = int(time / (3600 * 24))
hours = int((time % (3600 * 24)) / 3600)
minutes = int((time % 3600) / 60)
seconds = time % 60

print("time: %d days, %d hours, %d minutes & %d seconds" % (days, hours, minutes, seconds))
print("script made %.1f million calculations" % round((time / (delta_t * (10 ** 6))), 1))
print("script took %d seconds" % (t.time() - start))
print("--- Task A")
#1
a = 2
b = 4
c = a + b

values = [a, b, c]
#for x in values:
    #print(x)

a = a + 1
b = a
c = a + b

values = [a, b, c]
#for x in values:
    #print(x)

#2
x, y = (-11, -3)

z = (3*x) + (y**2)
print(z)

print("--- Task B")
#1
var = 3.14
varcopy = var

#2
MyPints = 3
drink2more = 2

MyPints += drink2more
MyPints += drink2more

#3
Num, Den = (3, 4)

Res = Num / Den

MyPints -= drink2more
print(MyPints)

print("--- Task C")
#1
a, b = (10, 5)
b, a = (a, b)

#2
c = 20

a, c = (c, a)
b, c = (c, b)

print(a)

print("--- Task D")
#1
b = 2
a = str(b)

#2
f = [a, b]

#3
c = "3"
d = int(c)

#4
g = [a, c]
h = [b, d]

#5
m = int(a) + int(c)

print(f)

from matplotlib import pyplot as plot

def arr_max(arr):
    m = None
    for x in arr:
        if x > m or m == None:
            m = x
    return m

def arr_min(arr):
    m = None
    for x in arr:
        if x < m or m == None:
            m = x
    return m

def arr_sum(arr):
    s = 0
    for x in arr:
        if x > s or s == None:
            s += x
    return s

print("--- Task A")

def trace(m):
    s = 0
    for i in range(len(m)):
        s += m[i][i]
    return s

print("--- Task B")

def mat_product(a, b):
    if len(a[0]) != len(b):
        return None
    n, m = len(a), len(b[0])
    c = []
    for j in range(m):
        row = []
        for i in range(n):
            row.append(arr_sum(list(
                map(lambda x, y: x*y[i], a[j], b)
            )))
        c.append(row)
    return c

print("--- Task C")

def is_anagram(a, b):
    for char in a:
        if char in b:
            b.remove(char)
        else:
            return False
    if len(b) == 0:
        return True
    else:
        return False

print("--- Task D")

with open("./exercise9/Temperatures.txt", "r") as f:
    t = [int(x) for x in f.read().split("\n")]

def temperature_analysis(t):
    i = 0
    t_avg = []
    max_t, min_t = None, None
    while(i < len(t)):
        t_daily = t[i:i+24]
        if arr_max(t_daily) > max_t or max_t == None:
            max_t = arr_max(t_daily)
        elif arr_min(t_daily) < min_t or min_t == None:
            min_t = arr_min(t_daily)
        t_avg.append(arr_sum(t_daily)/24)
        i += 24
    return t_avg, max_t, min_t

analysed_t, max_t, min_t = temperature_analysis(t)

plot.plot(range(len(analysed_t)), analysed_t)

plot.xlabel("Days")
plot.ylabel("Average Temperature")

plot.show()
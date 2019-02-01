from matplotlib import pyplot as plot
from re import compile, findall

print("--- Task A")

with open('./exercise4/CIDs.txt', 'r') as f:
    cid_text = f.read()

with open('./exercise4/Marks.txt', 'r') as f:
    marks_text = f.read()

cid_regex = marks_regex = compile('(\d+)\s+')

marks_list = [int(m) for m in marks_regex.findall(marks_text)]
cids_list = [int(cid) for cid in cid_regex.findall(cid_text)]

if len(marks_list) != len(cids_list):
    raise Exception("Lengths of cids and marks are not equal")

sum_marks = 0
largest_mark = 0

for m in marks_list[:60]:
    sum_marks += m
    if m > largest_mark:
        largest_mark = m

average_mark = sum_marks / len(marks_list[:60])
print("Average mark: %f, Largest mark: %i" % (average_mark, largest_mark))

print("--- Task B")

def find_score(cid):
    return marks_list[cids_list.index(cid)]

best_text = open('./exercise4/Best.txt', 'w+')

def return_cids(score):
    i = 0
    cids = []
    while i < len(marks_list):
        try:
            i += marks_list[i+1:].index(score) + 1
        except:
            break
        cids.append(str(cids_list[i]))
    return cids

best_text.write("\n".join(return_cids(largest_mark)))
best_text.close()

print(find_score(4597))

print("--- Task C")

def series_expansion(domain, n, step):
    x = domain[0]
    data_points = {'x':[], 'y':[]}
    while x <= domain[1]:
        y = 0
        for p in range(n+1):
            y += x ** p

        data_points['x'].append(x)
        data_points['y'].append(y)
        
        x += step
    
    return data_points

series_n2 = series_expansion([-0.8, 0.8], 2, 0.01)
series_n6 = series_expansion([-0.8, 0.8], 6, 0.01)
series_n10 = series_expansion([-0.8, 0.8], 10, 0.01)
series_n14 = series_expansion([-0.8, 0.8], 14, 0.01)

def reciprocal_one_minus_x():
    x = -0.8
    data_points = {'x':[], 'y':[]}
    while x <= 0.8:
        y = 1 / (1 - x)

        data_points['x'].append(x)
        data_points['y'].append(y)
        x += 0.01

    return data_points

y_x = reciprocal_one_minus_x()

plot.plot(
    series_n2['x'], series_n2['y'],
    series_n6['x'], series_n6['y'],
    series_n10['x'], series_n10['y'],
    series_n14['x'], series_n14['y'],
    y_x['x'], y_x['y'], 'r'
)

plot.show()

print(series_n14['y'][series_n14['x'].index(-0.70)])

print("--- Task D")

def approximate_reciprocal_one_minus_x(x, q):
    accurate = False
    n = 0
    y = 1
    while not accurate:
        n += 1
        new_y = y + (x ** n)

        if abs(new_y - y) < 0.1 ** q:
            accurate = True
        else:
            y = new_y

    true_y = 1 / (1 - x)
    print("Approximation: %f, Error: %f for %i elements in the series" % (y, abs(true_y - y), n+1))

approximate_reciprocal_one_minus_x(0.5, 6)
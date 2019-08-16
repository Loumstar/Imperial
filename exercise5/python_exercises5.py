from matplotlib import pyplot as plot

with open('./Groups.txt') as f:
    groups = f.read().split('\n')

with open('./Marks.txt') as f:
    marks = [int(a) for a in f.read().split('\n')]

with open('./Names.txt') as f:
    names = f.read().split('\n')

if not len(groups) == len(marks) == len(names):
    raise Exception("Length of text files are not equal")

print("--- Task A")

mecheng_tuple = [(x, groups[i], marks[i]) for i,x in enumerate(names)]

print(mecheng_tuple[10])

print("--- Task B")

def get_mark(s):
    return s[2]

def own_sort(l):
    sorted_list = [l[0]]
    for x in l[1:]:
        smallest = True
        for i,y in enumerate(sorted_list):
            if x[2] > y[2]:
                smallest = False
                sorted_list.insert(i, x)
                break
        if smallest:
            sorted_list.append(x)
    
    return sorted_list

mecheng_tuple_by_mark = sorted(mecheng_tuple, reverse=True, key=get_mark)


print("--- Task C")

mark_dict = {}
for student in mecheng_tuple:
    m = student[2]
    if m in mark_dict.keys():
        mark_dict[m][0] += 1
        mark_dict[m][1].append(student[0])
    else:
        mark_dict[m] = [1, [student[0]]]

mark_tuple = [(mark, info[0], info[1]) for mark,info in mark_dict.items()]

mark_vals = [m[0] for m in mark_tuple]
mark_counts = [m[1] for m in mark_tuple]

plot.plot(mark_vals, mark_counts, '*')

plot.xlabel("Mark Values")
plot.ylabel("Mark Counts")

plot.show()
print("--- Task A")

def read_matrix(lst, n, m):
    matrix = []
    for j in range(m):
        i = j * n
        matrix.append(lst[i:i+n])
    return matrix

def matrix_vector_dot(m, v):
    c = []
    for r in m:
        if len(r) != len(v):
            return 0
        c.append(sum(list(map(lambda x, y: x * y, r, v))))
    return c

with open("./exercise8/Q1A.txt", "r") as f:
    q1a_list = [int(x) for x in f.read().split("\n")]

q1a_matrix = read_matrix(q1a_list, 400, 200)

with open("./exercise8/Q1B.txt", "r") as f:
    q1b_list = [int(x) for x in f.read().split("\n")]
    
q1b_vector = read_matrix(q1b_list, 400, 1)[0]

print(matrix_vector_dot(q1a_matrix, q1b_vector)[34])

print("--- Task B")

def chess_board(n):
    board = []
    for y in range(n):
        row = [0, 1] * (n-1)
        row.insert(0, 1) if (y % 2 == 0) else row.append(0)
        board.append(row)
    return board

def write_chess_board_to_txt(n):
    f = open('./exercise8/chess_board.txt', 'w')
    f.write("11\n11\n")
    m = chess_board(n)
    for i in m:
        for ij in i:
            f.write("%d\n" % ij)
    f.close()

write_chess_board_to_txt(11)

print("--- Task C")

with open('./exercise8/Game1.txt', 'r') as f:
    list_g1 = [int(x) for x in f.read().split("\n")][2:]

g1 = read_matrix(list_g1, 8, 8)

with open('./exercise8/Game2.txt', 'r') as f:
    list_g2 = [int(x) for x in f.read().split("\n")][2:]

g2 = read_matrix(list_g2, 8, 8)

with open('./exercise8/Game3.txt', 'r') as f:
    list_g3 = [int(x) for x in f.read().split("\n")][2:]

g3 = read_matrix(list_g3, 8, 8)

with open('./exercise8/ChessA.txt', 'r') as f:
    list_q2 = [int(x) for x in f.read().split("\n")][2:]

q2 = read_matrix(list_q2, 8, 8)

def find_piece(m, val):
    for i in range(len(m)):
        for j in range(len(m[i])):
            if m[i][j] == val:
                return (i, j)
    return None

def check_adjacent_cells(m, val):
    c = find_piece(m, val)
    s = 0
    enemies = []
    for row in range(c[0], c[0]+2):
            if row in (0, len(m)):
                continue
            for col in range(c[1], c[1]+2):
                if col in (0, len(m[row])):
                    continue
                if m[row][col] * m[c[0]][c[1]] < 0:
                    print("Enemy adjacent to cell: (%d, %d) is %d" % (row, col, m[row][col]))
                    s += 1
                    enemies.append((row, col))
    print("There are %d enemies adjacent to the cell" % s)
    return enemies


def print_matrix(m):
    for r in m:
        print(r)

print(find_piece(q2, 1))
print(check_adjacent_cells(q2, 1))

print("--- Task D")

def transpose(m):
    r = []
    for j in range(len(m[0])):
        r.append(list(map(lambda x: x[j], m)))
    return r

with open('./exercise8/Trans.txt', 'r') as f:
    t_list = [int(x) for x in f.read().split("\n")][2:]

t = transpose(read_matrix(t_list, 402, 245))
print(t[11][25])
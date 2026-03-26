from collections import deque

import numpy as np
from scipy import linalg


def read_graph(filename):
    with open(filename, "r") as f:
        temp = f.readline()
        if not temp:
            return
        n = int(temp.strip())
        lst = [[] for i in range(n)]
        for i in range(n):
            line = f.readline().strip()
            if line:
                lst[i] = list(map(int, line.split()))
            else:
                lst[i] = []
    return lst


def bfs_distances(lst):  # расстояние от i до всех остальных вершин
    n = len(lst)
    paths = []
    for i in range(n):
        paths.append([])
        for j in range(n):
            paths[i].append([])
    distance_count = []
    bfs_order = []
    dist = [[0] * n for _ in range(n)]
    for i in range(n):
        visited = [False] * n
        q = deque([i])
        visited[i] = True
        order = []
        cnt = [0] * (n - 1)
        while q:
            u = q.popleft()
            order.append(u)
            for v in lst[u]:
                if not visited[v]:
                    dist[i][v] = dist[i][u] + 1
                    dist[v][i] = dist[i][v]
                    cnt[dist[i][v] - 1] += 1
                    visited[v] = True
                    q.append(v)
                    paths[i][v] = paths[i][u].copy() + [v]
        bfs_order.append(order)
        distance_count.append(cnt)
    return dist, paths, bfs_order, distance_count


def get_binary(n, p):
    if len(p) == 0:
        return np.array([1])
    N = 2**n
    b = (np.arange(N)[:, None] >> np.arange(n - 1, -1, -1)) & 1
    # ТО ЖЕ ЧТО И 1 - 2*(np.sum(b[:,p], axis=1) % 2)
    return 1 - (
        (np.sum(b[:, p], axis=1) & 1) << 1
    )  # возвращаем список единиц и минус единиц


def get_I(k):
    return np.eye(2 ** (k))


# lst = read_graph("falcon_r5_11h.txt")
lst = read_graph("graph1.txt")
n = len(lst)
dist, paths, bfs_order, distance_count = bfs_distances(lst)
# paths без начальной вершины
M_CNOT = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
X = np.array([[0, 1], [1, 0]])
k = n - 1
V_LIST = np.zeros((n - 1, 2 ** (n - 1), 2 ** (n - 1)))
print()
V_LIST[0] = np.eye(2 ** (n - 1))  # массив V
for k in range(1, n - 1):
    I_first = get_I(n - k - 1)
    I_second = get_I(n - k - 2)
    X_repeats = np.zeros((2 ** (n - k - 1), 2 ** (n - k - 1)))
    X_repeats[2 ** (n - k - 2) :, 0 : 2 ** (n - k - 2)] = I_second
    X_repeats[0 : 2 ** (n - k - 2), 2 ** (n - k - 2) :] = I_second
    temp = np.concatenate((I_first[None, :, :], X_repeats[None, :, :]), axis=0)
    V_LIST[k] = linalg.block_diag(
        *np.tile(temp, (2 ** (k - 1), 1, 1))
    )  # здесь был test

    # np.tile (A, a,b,c) A - список матриц, a,b,c соответствуют повторениям на других осях (сколько осей столько и букв!)

    # V_LIST[k]= np.kron(np.kron(get_I(k - 1), M_CNOT), get_I(n - k - 2)) # подозреваю, что это правильно, но медленно (САМАЯ ВЕРХНЯЯ ФОРМУЛА НА 9-ОЙ СТРАНИЦЕ)
    # assert(np.all(V_LIST[k] == test))

    print(V_LIST[k], end="\n\n")

print("=========================================")
H_volna = np.array([[1, 1], [1, -1]])
U_LIST = np.zeros((n - 1, 2 ** (n - 1), 2 ** (n - 1)))
for k in range(n - 1):
    I = get_I(n - k - 2)
    I_repeats = np.tile(I, (2, 2))
    I_repeats[2 ** (n - k - 2) :, 2 ** (n - k - 2) :] = -I_repeats[
        2 ** (n - k - 2) :, 2 ** (n - k - 2) :
    ]
    U_LIST[k] = 0.5 * linalg.block_diag(
        *np.tile(I_repeats[None, :, :], (2**k, 1, 1))
    )  # здесь был test
    pass

    # U_LIST[k] = 0.5 * np.kron(np.kron(get_I(k), H_volna), get_I(n - k - 2)) # не оптимальный вариант
    # assert(np.all(U_LIST[k] == test))

    # print(U_LIST[k], end='\n\n')

Q_LIST = np.zeros((n - 1, 2 ** (n - 1), 2 ** (n - 1)))

for k in range(n - 1):
    p = np.array(paths[k][n - 1][:-1]) - (
        k + 1
    )  # преобразованный путь между таргетом и контролирующим битом на k-ом шаге
    d = get_binary(n - k - 2, p)
    Q_LIST[k] = np.diag(np.repeat(np.concatenate((np.ones(2 ** (n - k - 2)), d)), 2**k))

V_multiplies = np.linalg.multi_dot(V_LIST[::-1])
print(V_multiplies)
MULT_QU = Q_LIST @ U_LIST
MULT_QU = np.linalg.multi_dot(MULT_QU[::-1])
W_LIST = V_multiplies @ MULT_QU
print(W_LIST)

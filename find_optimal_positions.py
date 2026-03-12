import time
from collections import deque

import matplotlib.pyplot as plt
from qiskit import QuantumCircuit
from qiskit.primitives import StatevectorSampler
from qiskit.visualization import circuit_drawer, plot_histogram


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


def floyd_warshall(lst):
    n = len(lst)
    dist = [[11111111] * n for _ in range(n)]
    for i in range(n):
        dist[i][i] = 0
        for v in lst[i]:
            dist[i][v] = 1
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if dist[i][k] + dist[k][j] < dist[i][j]:
                    dist[i][j] = dist[i][k] + dist[k][j]
    return dist


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


def gray_changed_bits_cyclic(n):
    gc = [i ^ (i >> 1) for i in range(1 << n)]
    changed = [(gc[i] ^ gc[(i + 1) % (1 << n)]).bit_length() - 1 for i in range(1 << n)]
    return changed


def gray_changed_bits_cyclic_reverse(n):
    """
    Возвращает список индексов битов (слева направо), которые меняются при переходе
    от одного числа к следующему в Gray code длины n, циклически.
    """
    gc = [i ^ (i >> 1) for i in range(1 << n)]
    changed = []
    for i in range(1 << n):
        diff = gc[i] ^ gc[(i + 1) % (1 << n)]
        # номер бита слева направо, нумерация С ЕДИНИЦЫ
        bit_index = n - (diff.bit_length() - 1) - 1
        changed.append(bit_index)
    return changed


# поиск таргета, с учетом расположения всех остальных вершин
def get_cnot_with_optimal_location(distance_count, n):
    lst = [0] * n
    for i in range(n):
        cnt = 0
        for j in range(n - 1):
            for t in range(distance_count[i][j]):
                if (
                    cnt < n - 2
                ):  # два последних слоя кодов Грея встречаются одинаково часто
                    cnt += 1
                    lst[i] += 2 ** (n - 1 - (t + 1 + sum(distance_count[i][:j]))) * (
                        2 * (j + 1) - 1
                    )
                else:
                    lst[i] += 2 * (2 * (j + 1) - 1)
    # print(*lst, sep='\n\n')
    return lst


def get_cnot_not_optimal_location(gray_codes, dist, n):
    # поиск таргетов, для которых сумма расстояний от всех вершин до таргета минимальна. (необходимое условие, но не достаточное)
    number_cnots = [0] * n
    r = len(gray_codes)
    for target in range(n):
        for i in range(r // 2):
            control_bit = gray_codes[i]
            if control_bit == target:
                control_bit = n - 1
            temp = 2 * (dist[target][control_bit]) - 1
            number_cnots[target] += temp
        number_cnots[target] = 2 * number_cnots[target]
    return number_cnots


def draw_optimal_circuit(good_target, gray_codes, optimal_order, paths, n):
    qc = QuantumCircuit(n)
    graycode_vertex = (
        {}
    )  # словарь, где ключом является код изменившегося бита в таблице кодов Грея,
    # а значение - соответствующая ему вершина в графе.
    r = len(gray_codes)
    for k in range(r // 2):
        temp = gray_codes[k]
        if temp not in graycode_vertex:
            graycode_vertex[temp] = optimal_order[-temp - 1]
    print(graycode_vertex)
    print("optimal_order", optimal_order)
    vertex_qubit = {}
    for i in range(len(optimal_order)):
        vertex_qubit[optimal_order[i]] = n - 1 - i
    print(vertex_qubit)
    for k in range(r):
        qc.ry(0, n - 1)  # на последнем месте в квантовой схеме всегда target
        qc.barrier()
        control_bit = graycode_vertex[
            gray_codes[k]
        ]  # вершина графа, соответствующая коду из таблицы Грея
        qc.cx(gray_codes[k], vertex_qubit[paths[control_bit][good_target][0]])
        for j in range(len(paths[control_bit][good_target]) - 1):
            qc.cx(
                vertex_qubit[paths[control_bit][good_target][j]],
                vertex_qubit[paths[control_bit][good_target][j + 1]],
            )
        for j in range(len(paths[good_target][control_bit]) - 1):
            qc.cx(
                vertex_qubit[paths[good_target][control_bit][j + 1]],
                vertex_qubit[paths[good_target][control_bit][j]],
            )
        qc.barrier()
        # print(circuit_drawer(qc, output="text"))

    # print(bfs_order[good_target])

    circuit_drawer(qc, output="mpl")
    plt.show()


def draw_not_optimal_circuit(
    good_target, gray_codes, paths, n
):  # отрисовка схемы, где мы меняем между собой target и один из контролирующих
    r = len(gray_codes)
    qc = QuantumCircuit(n)
    for k in range(r):
        control_bit = gray_codes[k]
        if control_bit == good_target:
            control_bit = n - 1
        qc.ry(0, good_target)
        qc.barrier()
        qc.cx(control_bit, paths[control_bit][good_target][0])
        for j in range(len(paths[control_bit][good_target]) - 1):
            qc.cx(
                paths[control_bit][good_target][j],
                paths[control_bit][good_target][j + 1],
            )
        for j in range(len(paths[good_target][control_bit]) - 1):
            qc.cx(
                paths[good_target][control_bit][j + 1],
                paths[good_target][control_bit][j],
            )
        qc.barrier()
        print(circuit_drawer(qc, output="text"))
    circuit_drawer(qc, output="mpl")
    plt.show()


lst = read_graph("graph1.txt")
# print(lst)
# print(floyd_warshall(lst))
# print(bfs_distances(lst))

# сравнение Флойда и BFS
# start = time.perf_counter()
# fw_dist = floyd_warshall(lst)
# end = time.perf_counter()
# print("Floyd-Warshall:", end - start, "секунд")

# start = time.perf_counter()
# bfs_dist = bfs_distances(lst)
# end = time.perf_counter()
# print("BFS:", end - start, "секунд")

dist, paths, bfs_order, distance_count = bfs_distances(lst)
n = len(lst)
optimal_cnots = get_cnot_with_optimal_location(distance_count, n)
# print(optimal_cnots)
gray_codes = gray_changed_bits_cyclic_reverse(n - 1)

number_cnots = get_cnot_not_optimal_location(gray_codes, dist, n)
good_not_optimal_target = number_cnots.index(min(number_cnots))

# draw_not_optimal_circuit(good_not_optimal_target, gray_codes, paths, n)

good_optimal_target = optimal_cnots.index(
    min(optimal_cnots)
)  # номер вершины графа, которую можно взять за target
optimal_order = bfs_order[good_optimal_target]
draw_optimal_circuit(good_optimal_target, gray_codes, optimal_order, paths, n)

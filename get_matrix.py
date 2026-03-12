import numpy as np


def get_matrix_b_g(n):
    N = 2**n
    b = (np.arange(N)[:, None] >> np.arange(n - 1, -1, -1)) & 1
    g = (
        (np.arange(N)[:, None] ^ (np.arange(N)[:, None] >> 1))
        >> np.arange(n - 1, -1, -1)
    ) & 1
    t = (b @ g.T) % 2
    return (-1) ** t


# n - количество контролирующих кубитов, n + 1 всего кубитов
n = int(input())
teta = np.random.uniform(0, 2 * np.pi, size=(2**n, 1))
matrix_M = get_matrix_b_g(n)
print(matrix_M)
alpha = np.pow(-n, 2) * (matrix_M.T @ teta)
print(alpha, end="\n\n")
alpha_mod = np.mod(alpha, 2 * np.pi)
print(alpha_mod)

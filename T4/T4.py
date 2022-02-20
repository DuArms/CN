import time

from T3.T3 import *

EPS = 1e-14


def tridiagonal_main_diagonal_not_null(tridiagonal):
    for x in tridiagonal[0]:
        if abs(x) < EPS:
            return False
    return True


def read_free_terms(file_name="T4/f1.txt"):
    with open(file_name, "r") as file:
        lines = file.readlines()

        free_terms = []
        n = int(lines[0])

        for i in range(2, n + 2):
            free_terms.append(float(lines[i]))

        return free_terms


def tridiagonal_solve_ecuation(A, F):
    n = len(F)
    p, q = n - len(A[1]), n - len(A[2])

    X = [2 for x in range(n)]

    start_time = time.time()
    while True:
        X_sum_old = sum(X)

        for i in range(n):
            x = F[i]

            # for j in {i - p, i, i + q}:
            #     if 0 <= j < i:
            #         x -= tridiagonal_matrix_get_el(A, i, j) * X[j]

            if 0 <= i - p < i:
                x -= tridiagonal_matrix_get_el(A, i, i - p) * X[i - p]

            # for j in {i - p, i, i + q}:
            #     if i < j < n:
            #         x -= tridiagonal_matrix_get_el(A, i, j) * X[j]

            if i < i + q < n:
                x -= tridiagonal_matrix_get_el(A, i, i + q) * X[i + q]

            X[i] = x / tridiagonal_matrix_get_el(A, i, i)

        if abs(sum(X) - X_sum_old) <= EPS:
            break

        if time.time() - start_time > 1:
            break

    return X


def tridiagonal_solution_norm(A, X, F):
    n = len(F)
    p, q = n - len(A[1]), n - len(A[2])
    norm = None

    for i in range(n):

        s = 0
        for j in {i - p, i, i + q}:
            if 0 <= j < n:
                s += tridiagonal_matrix_get_el(A, i, j) * X[j]

        s -= F[i]
        s = abs(s)

        if norm is None:
            norm = s
        else:
            norm = max(norm, s)

    return norm


def homework_4_presentation():
    nr_tests = 6

    tridiagonals = [[] for x in range(nr_tests)]
    for i in range(nr_tests):
        tridiagonals[i] = read_triangular_matrix("T4/a" + str(i+1) + ".txt")

    for i in range(nr_tests):
        print(f"Tridiagonal[{i+1}]'s main diagonal: {tridiagonal_main_diagonal_not_null(tridiagonals[i])}")

    free_terms = [[] for x in range(nr_tests)]
    for i in range(nr_tests):
        free_terms[i] = read_free_terms("T4/f" + str(i+1) + ".txt")

    solutions = [[] for x in range(nr_tests)]
    for i in range(nr_tests):
        solutions[i] = tridiagonal_solve_ecuation(tridiagonals[i], free_terms[i])

    for i in range(nr_tests):
        print(f"Ecuation {i+1} norm: ", end="")
        print(tridiagonal_solution_norm(tridiagonals[i], solutions[i], free_terms[i]))



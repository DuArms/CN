import numpy as np
import copy
import sys
from scipy.linalg import lu

EPS = 1e-12


def generate_spd_matrix(N=3):
    matrix = np.array([np.random.uniform(0, 1, N) for _ in range(1, N + 1)])
    matrix = matrix @ matrix + np.eye(N) * N
    return matrix, np.random.sample(N)


def read_matrix_from_cmd():
    N = int(input("N = "))
    print("Enter A matrix :")
    matrix = [(list(map(float, input(f" Line[{_}] : ").split()))) for _ in range(N)]
    print(f"Enter B matrix ({N},1) :")
    b_matrix = [(list(map(float, input(f" Line[{_}] : ").split()))[0]) for _ in range(N)]
    return matrix, b_matrix


def read_matrix_from_file(input_name="T2/T2_1.in"):
    with open(input_name, "r") as file:
        lines = file.readlines()
        n = int(lines[0])

        A = []
        for i in range(n):
            row = []
            for w in lines[1 + i].split():
                row.append(float(w))
            A.append(row)

        B = []
        for i in range(n):
            B.append(float(lines[1 + n + i]))

        return A, B


def cholesky_factorization(A):
    return cholesky_factorization_bonus(A)

    L = []
    n = len(A)
    for p in range(n):
        row = [0 for p in range(len(A[p]))]
        L.append(row)

    for p in range(len(A)):
        s = 0
        for j in range(p):
            s += L[p][j] ** 2

        if A[p][p] - s < -EPS:
            return False, -1

        L[p][p] = np.sqrt(A[p][p] - s)

        for i in range(p + 1, n):
            m = 0
            for j in range(p):
                m += L[i][j] * L[p][j]
            L[i][p] = (A[i][p] - m) / L[p][p]

    return True, L


def cholesky_factorization_bonus(A_init):
    A = copy.deepcopy(A_init)
    n = len(A)

    d = [A[i][i] for i in range(n)]

    for p in range(len(A)):
        s = 0
        for j in range(p):
            s += A[p][j] ** 2

        if d[p] - s < -EPS:
            return False, -1

        A[p][p] = np.sqrt(d[p] - s)

        for i in range(p + 1, n):
            m = 0
            for j in range(p):
                m += A[i][j] * A[p][j]
            A[i][p] = (A[i][p] - m) / A[p][p]

    return True, A


def substitution_method_down(L, B):
    n = len(L)

    Ys = []
    for i in range(n):
        s = 0
        for j in range(i):
            s += Ys[j] * L[i][j]
        Ys.append((B[i] - s) / L[i][i])

    return Ys


def substitution_method_up(L, B):
    n = len(L)

    Xs = [0 for i in range(n)]
    for i in range(n - 1, -1, -1):
        s = 0
        for j in range(i + 1, n):
            s += L[j][i] * Xs[j]
        Xs[i] = (B[i] - s) / L[i][i]

    return Xs


def check_solution(A, B, X):
    n = len(A)

    for y in range(n):
        s = 0

        for i in range(n):
            if i <= y:
                Z = A[y][i]
            else:
                Z = A[i][y]
            s += Z * X[i]

        print("\nLine[" + str(y) + "]: ")
        print("\t\tCalculations made with our solution:\t" + str(s))
        print("\t\tExpected Result:\t" + str(B[y]))


def cholesky_determinant(L):
    n = len(L)
    prod = 1
    for i in range(n):
        prod *= L[i][i]
    return prod


def convert_to_full_matrix_duplicate(L):
    n = len(L)
    A = [[0 for i in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if j <= i:
                Z = L[i][j]
            else:
                Z = L[j][i]
            A[i][j] = Z
    return A


def convert_to_full_matrix_half(L):
    n = len(L)
    A = [[0 for i in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            if j <= i:
                Z = L[i][j]
            else:
                continue
            A[i][j] = Z
    return A


def matrix_norm(A, X, B):
    n = len(A)
    norm = 0
    for y in range(n):
        s = 0

        for i in range(n):
            if i <= y:
                Z = A[y][i]
            else:
                Z = A[i][y]
            s += Z * X[i]

        norm += (s - B[y]) ** 2

    return norm ** (1 / 2)


def lower_upper_decomposition(A):
    n = len(A)
    U = [[0 for i in range(n)] for j in range(n)]
    L = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        for k in range(i, n):
            s = 0
            for j in range(i):
                s += L[i][j] * U[j][k]

            U[i][k] = A[i][k] - s

        for k in range(i, n):
            if i == k:
                L[k][i] = 1
            else:
                s = 0
                for j in range(i):
                    s += L[k][j] * U[j][i]

                L[k][i] = (A[k][i] - s) / U[i][i]

    return L, U


def print_tm(A):
    if type(A[0]) not in {list, np.ndarray}:
        A = np.array(A).reshape((len(A), 1))
    A = np.array(A)
    print(A)


def print_named_matrix(name, matrix):
    print(f"{name}:")
    print_tm(matrix)
    print()


def homework_2_presentation():
    # A, B = read_matrix_from_cmd()
    A, B = read_matrix_from_file()
    # A, B = generate_spd_matrix()

    A_full = convert_to_full_matrix_duplicate(A)
    n = len(A)

    print("Matrix entered A : ")
    print_tm(A_full)
    print("Matrix entered B : ")
    print_tm(B)

    possible, L = cholesky_factorization(A)

    if not possible:
        print("Matrix A is not positive defined")
        return

    print("\nCholesky factorization :")
    Lc = convert_to_full_matrix_half(L)
    Lt = [[Lc[i][j] for i in range(n)] for j in range(n)]
    print_named_matrix("L", Lc)
    print_named_matrix("Lt", Lt)

    Det_L = cholesky_determinant(L)
    print(f"\nDeterminant(L) Handmade: {Det_L}")
    print(f"Determinant(A) Handmade: {Det_L * Det_L}")

    print(f"\nDeterminant(L) NumPy: {np.linalg.det(convert_to_full_matrix_half(L))}")
    print(f"Determinant(A) NumPy: {np.linalg.det(A_full)}")

    Y = substitution_method_down(L, B)
    X = substitution_method_up(L, Y)

    X_lib = np.linalg.solve(A_full, B)

    print_named_matrix("\nSolution for LY = b", Y)
    print_named_matrix("Solution for LtX = Y", X)
    print_named_matrix("Solution from library", X_lib)

    print("\nCheck the validity of the solution:", end="")
    check_solution(A=A, B=B, X=X)

    norm = matrix_norm(A=A, X=X, B=B)
    print(f"\n||AX-B|| = {norm}\n")
    if abs(norm) < EPS:
        print("Found correct solutions!")
    else:
        print("Wrong solution!")

    L1, U1 = lower_upper_decomposition(A_full)
    p, L2, U2 = lu(A_full)

    print("Our LU:")
    print_named_matrix("\tL", L1)
    print_named_matrix("\tU", U1)
    print("\nSkiPy LU:")
    print_named_matrix("\tL", L2)
    print_named_matrix("\tU", U2)

    inv_a_chol = []
    for i in range(0,n):
        b = [0 for i in range(n)]
        b[i] = 1
        column = substitution_method_up(L, substitution_method_down(L, b))
        inv_a_chol.append(column)

    inv_a_chol = np.array(inv_a_chol).T

    inv_lt = np.linalg.inv(Lt)
    inv_a_google = inv_lt @ inv_lt.T
    inv_a_np = np.linalg.inv(A_full)

    print_named_matrix("\nInvers of A with numpy :", inv_a_np)
    print_named_matrix("Invers of A with google method:", inv_a_google)
    print_named_matrix("Invers of A with chol method:", inv_a_chol)
    norm = np.linalg.norm(inv_a_chol - inv_a_np)
    print(f"|| A_chol^-1  - A_lib^-1 || = {norm}")
    norm = np.linalg.norm(inv_a_google - inv_a_np)
    print(f"|| A_google^-1  - A_lib^-1 || = {norm}")

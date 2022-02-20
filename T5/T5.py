from T2.T2 import *
import numpy as np

EPS = 1e-6


def matrix_symetric_gel_el(A, i, j):
    if j >= i:
        return A[j][i]
    return A[i][j]


def v(i, j):
    if j <= i:
        return int(i*(i+1)/2 + j)
    return int(j*(j+1)/2 + i)


def iacobi_aprox(simetric_matrix, k_max=10000):
    A = simetric_matrix
    n = len(A)

    U = []
    for i in range(n):
        row = [0 for p in range(n)]
        row[i] = 1
        U.append(row)

    k = 0
    while True:

        if k >= k_max:
            break

        p, q = -1, -1
        max_v = -1

        for i in range(1, n):
            for j in range(i):
                if abs(A[i][j]) > EPS and abs(A[i][j]) > max_v:
                    p, q = i, j
                    max_v = abs(A[i][j])

        if max_v == -1:
            break

        A_ = lambda i, j: matrix_symetric_gel_el(A, i, j)

        alpha = (A_(p, p) - A_(q, q)) / (2 * A_(p, q))
        if alpha >= 0:
            t = -alpha + np.sqrt(alpha * alpha + 1)
        else:
            t = -alpha - np.sqrt(alpha * alpha + 1)

        c = 1.0 / np.sqrt(1 + t * t)
        s = t / np.sqrt(1 + t * t)

        B = []
        for i in range(n):
            row = [A[i][p] for p in range(len(A[i]))]
            B.append(row)

        B[p][p] = c * c * A_(p, p) + s * s * A_(q, q) + 2 * c * s * A_(p, q)
        B[q][q] = s * s * A_(p, p) + c * c * A_(q, q) - 2 * c * s * A_(p, q)
        if q <= p:
            B[p][q] = (c * c - s * s) * A_(p, q) + c * s * (A_(q, q) - A_(p, p))
        else:
            B[q][p] = (c * c - s * s) * A_(p, q) + c * s * (A_(q, q) - A_(p, p))
        for j in range(n):
            if j != p and j != q and j <= p:
                B[p][j] = c * A_(p, j) + s * A_(q, j)

            if j != p and j != q and p <= j:
                B[j][p] = c * A_(p, j) + s * A_(q, j)

            if j != p and j != q and j <= q:
                B[q][j] = -s * A_(p, j) + c * A_(q, j)

            if j != p and j != q and q <= j:
                B[j][q] = -s * A_(p, j) + c * A_(q, j)
        A = B

        V = []
        for i in range(n):
            row = [U[i][p] for p in range(n)]
            V.append(row)

        for i in range(n):
            V[i][p] = c * U[i][p] + s * U[i][q]
            V[i][q] = -s * U[i][p] + c * U[i][q]
        U = V

        k += 1

    return A, U


def iacobi_aprox_bonus(simetric_matrix, k_max=10000):
    A = simetric_matrix
    n = len(A)

    # Deconstruction of A[n, n] -> A[n*(n+1)/2]
    B = []
    for i in range(n):
        for j in range(i+1):
            B.append(A[i][j])
    A = B

    U = []
    for i in range(n):
        row = [0 for p in range(n)]
        row[i] = 1
        U.append(row)

    k = 0
    while True:

        if k >= k_max:
            break

        p, q = -1, -1
        max_v = -1

        for i in range(1, n):
            for j in range(i):
                if abs(A[v(i, j)]) > EPS and abs(A[v(i, j)]) > max_v:
                    p, q = i, j
                    max_v = abs(A[v(i, j)])

        if max_v == -1:
            break

        alpha = (A[v(p, p)] - A[v(q, q)]) / (2 * A[v(p, q)])
        if alpha >= 0:
            t = -alpha + np.sqrt(alpha * alpha + 1)
        else:
            t = -alpha - np.sqrt(alpha * alpha + 1)

        c = 1.0 / np.sqrt(1 + t * t)
        s = t / np.sqrt(1 + t * t)

        B = []
        for el in A:
            B.append(el)

        B[v(p, p)] = c * c * A[v(p, p)] + s * s * A[v(q, q)] + 2 * c * s * A[v(p, q)]
        B[v(q, q)] = s * s * A[v(p, p)] + c * c * A[v(q, q)] - 2 * c * s * A[v(p, q)]

        B[v(p, q)] = (c * c - s * s) * A[v(p, q)] + c * s * (A[v(q, q)] - A[v(p, p)])

        for j in range(n):
            if j != p and j != q:
                B[v(p, j)] = c * A[v(p, j)] + s * A[v(q, j)]

            if j != p and j != q:
                B[v(q, j)] = -s * A[v(p, j)] + c * A[v(q, j)]
        A = B

        V = []
        for i in range(n):
            row = [U[i][p] for p in range(n)]
            V.append(row)

        for i in range(n):
            V[i][p] = c * U[i][p] + s * U[i][q]
            V[i][q] = -s * U[i][p] + c * U[i][q]
        U = V

        k += 1

    B = []
    for i in range(n):
        row = []
        for j in range(i+1):
            row.append(A[v(i, j)])
        B.append(row)

    return B, U


def cholesky_decomposition_convergence(A, k_max=10000):
    A_new = A
    for k in range(k_max):
        Junk, L = cholesky_factorization(A)
        Lnp = np.array(convert_to_full_matrix_half(L))
        A_new = Lnp.T @ Lnp

        if np.linalg.norm(A-A_new) < EPS:
            break

        A = A_new

    return A_new


def generate_random_matrix(p, n):
    matrix = np.array([np.random.uniform(0, 1, n) for _ in range(0, p)])
    return matrix


def homework_5_presentation():

    A, junk = read_matrix_from_file()
    # A, junk = generate_spd_matrix()
    print_named_matrix("A", convert_to_full_matrix_duplicate(A))

    # Jacobi Convergence
    print("Jacobi Convergence")

    B, U = iacobi_aprox_bonus(A)
    print_named_matrix("U", U)

    print_named_matrix("A Final", convert_to_full_matrix_duplicate(B))

    Unp = np.array(U)
    Anp = np.array(convert_to_full_matrix_duplicate(A))
    print_named_matrix("U^T * A * U", Unp.T @ Anp @ Unp)

    print("||AU - UB||")
    Bnp = np.array(convert_to_full_matrix_duplicate(B))
    print(np.linalg.norm(Anp @ Unp - Unp @ Bnp))

    eigvals = [B[i][i] for i in range(len(A))]
    print("\nEigvals (Jacobi)")
    print(eigvals)

    print("\nEigvals (NumPy)")
    print(np.linalg.eigvals(convert_to_full_matrix_duplicate(A)))

    # Cholesky Convergence
    print("\n------------------------------------------------------------------")
    print("\nCholesky Convergence")
    A = convert_to_full_matrix_duplicate(A)
    print_named_matrix("A Init", A)

    A_K = cholesky_decomposition_convergence(A)
    print_named_matrix("A Final", A_K)

    # Singular Values
    print("\n------------------------------------------------------------------")
    p, n = 5, 4
    print("\nSingular Value Decomposition")
    B = generate_random_matrix(p=p, n=n)
    print_named_matrix("B", B)

    U, S, Vh = np.linalg.svd(B, full_matrices=True)
    print("Singular Values:", S)
    print("Rang(B): ", sum(map(lambda x: x > 0, S)))
    print("Conditioning Number(B): ", max(S) / min(filter(lambda x: x > 0, S)))

    SI = np.zeros(shape=(n, p))
    for i, x in enumerate(filter(lambda x: x > 0, S)):
        SI[i][i] = 1.0 / x

    BI = (Vh.T @ SI) @ U.T
    print_named_matrix("Moore Penrose Pseudo Inverse(B)", BI)

    BJ = np.linalg.inv(B.T @ B) @ B.T
    print_named_matrix("Smallest Squares Pseudo Inverse(B)", BJ)

    print("||BI - BJ|| = ", np.linalg.norm(BI-BJ))

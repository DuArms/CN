from bisect import bisect_left


def sparse_matrix_add_el(sparse_matrix, i, j, val):
    for x in sparse_matrix[i]:
        if x[0] == j:
            x[1] += val
    sparse_matrix[i].append([j, val])


def sparse_matrix_get_el(sparse_matrix, i, j):
    for x in sparse_matrix[i]:
        if x[0] == j:
            return x[1]
    return 0


def read_sparse_matrix(file_name="T3/a.txt"):
    with open(file_name, "r") as file:
        lines = file.readlines()
        n = int(lines[0])

        sparse_matrix = [[] for i in range(n)]
        for i in range(2, len(lines)):
            if lines[i] == "\n" or lines[i] == "":
                break

            lines[i] = lines[i].replace(",", "")
            words = lines[i].split()
            val = float(words[0])
            i = int(words[1])
            j = int(words[2])

            sparse_matrix_add_el(sparse_matrix, i, j, val)

        for i, line in enumerate(sparse_matrix):
            sparse_matrix[i] = sorted(line)

        while len(sparse_matrix[-1]) == 0:
            sparse_matrix.pop()

        return sparse_matrix


def read_triangular_matrix(file_name="T3/b.txt"):
    with open(file_name, "r") as file:
        lines = file.readlines()

        tridiagonal_matrix = [[], [], []]
        n = int(lines[0])
        p = int(lines[1])
        q = int(lines[2])

        for i in range(4, n + 4):
            tridiagonal_matrix[0].append(float(lines[i]))
        for i in range(n + 5, n + 5 + (n - p)):
            tridiagonal_matrix[1].append(float(lines[i]))
        for i in range(n + 6 + (n - p), n + 6 + (n - p) + (n - q)):
            tridiagonal_matrix[2].append(float(lines[i]))

        return tridiagonal_matrix


def addition(sparse_matrix, tridiagonal_matrix):
    n = len(sparse_matrix)
    p, q = n - len(tridiagonal_matrix[1]), n - len(tridiagonal_matrix[2])

    result = [[] for i in range(n)]
    for i in range(n):
        line = []

        added_p = False
        added_q = False
        added_n = False

        for j, el in sparse_matrix[i]:

            if i == j:
                line.append([j, el + tridiagonal_matrix[0][i]])
                added_n = True
            elif i - j == p:
                line.append([j, el + tridiagonal_matrix[1][i - p]])
                added_p = True
            elif j - i == q:
                line.append([j, el + tridiagonal_matrix[2][i]])
                added_q = True
            else:
                line.append([j, el])

        if not added_n:
            line.append([i, tridiagonal_matrix[0][i]])

        if not added_p and i >= p:
            line.append([i - p, tridiagonal_matrix[1][i - p]])

        if not added_q and i < n - q:
            line.append([i + q, tridiagonal_matrix[2][i]])

        result[i] = sorted(line)

    return result


def multiplication(sparse_matrix, tridiagonal_matrix):
    n = len(sparse_matrix)
    p, q = n - len(tridiagonal_matrix[1]), n - len(tridiagonal_matrix[2])

    result = [[] for i in range(n)]

    for i in range(n):
        to_hunt = set()
        for j, val in sparse_matrix[i]:
            if j - q >= 0:
                to_hunt.add(j - q)
            to_hunt.add(j)
            to_hunt.add(j + p)

        for j in to_hunt:
            s = 0
            if j < len(tridiagonal_matrix[0]):
                s = tridiagonal_matrix[0][j] * sparse_matrix_get_el(sparse_matrix, i, j)
            if j < len(tridiagonal_matrix[1]):
                s += tridiagonal_matrix[1][j] * sparse_matrix_get_el(sparse_matrix, i, j + p)
            if j - q < len(tridiagonal_matrix[2]):
                s += tridiagonal_matrix[2][j - q] * sparse_matrix_get_el(sparse_matrix, i, j - q)
            if s:
                sparse_matrix_add_el(result, i, j, s)

        result[i] = sorted(result[i])

    return result


def tridiagonal_matrix_get_el(tridiagonal, i, j):
    n = len(tridiagonal[0])
    p, q = n - len(tridiagonal[1]), n - len(tridiagonal[2])

    if i == j and i < n:
        return tridiagonal[0][i]
    if i-j == p and i >= p:
        return tridiagonal[1][i-p]
    if j-i == q and i < n - q:
        return tridiagonal[2][i]

    return 0


def multiplication_tridiagonals(tridiagonal1, tridiagonal2):
    n = len(tridiagonal1[0])
    p1, q1 = n - len(tridiagonal1[1]), n - len(tridiagonal1[2])
    p2, q2 = n - len(tridiagonal2[1]), n - len(tridiagonal2[2])

    result = [[] for i in range(n)]

    for i in range(n):
        to_hunt = set()
        for j in {i - p1, i, i + q1}:
            if n > j - q2 >= 0:
                to_hunt.add(j - q2)
            if n > j >= 0:
                to_hunt.add(j)
            if 0 <= j + p2 < n:
                to_hunt.add(j + p2)

        for j in to_hunt:
            s = 0

            if j < len(tridiagonal2[0]):
                s = tridiagonal2[0][j] * tridiagonal_matrix_get_el(tridiagonal1, i, j)
            if j < len(tridiagonal2[1]):
                s += tridiagonal2[1][j] * tridiagonal_matrix_get_el(tridiagonal1, i, j + p2)
            if j - q2 < len(tridiagonal2[2]):
                s += tridiagonal2[2][j - q2] * tridiagonal_matrix_get_el(tridiagonal1, i, j - q2)
            if s:
                sparse_matrix_add_el(result, i, j, s)

        result[i] = sorted(result[i])

    return result


def homework_3_presentation():
    sparse_matrix = read_sparse_matrix()
    tridiagonal_matrix = read_triangular_matrix()

    addition_res = addition(sparse_matrix, tridiagonal_matrix)
    addition_pre = read_sparse_matrix("T3/aplusb.txt")
    print(addition_res == addition_pre)

    multiplication_res = multiplication(sparse_matrix, tridiagonal_matrix)
    multiplication_pre = read_sparse_matrix("T3/aorib.txt")
    print(multiplication_res == multiplication_pre)

    tridiagonal1 = read_triangular_matrix("T3/c.txt")
    tridiagonal2 = read_triangular_matrix("T3/d.txt")

    multiplication_res = multiplication_tridiagonals(tridiagonal1, tridiagonal2)
    multiplication_pre = read_sparse_matrix("T3/corid.txt")
    print(multiplication_res == multiplication_pre)

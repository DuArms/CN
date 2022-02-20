import numpy as np
import time

PI = np.math.pi
PI_2 = np.math.pi / 2
PI_4 = np.math.pi / 4


def machine_precision_linear():
    x = 1.0
    u = 1.0

    while True:
        if x + u == x:
            break
        u /= 10

    return u


def machine_precision_binary():
    right = 1
    left = 0
    base = 10
    number = base

    while number + 1 != 1:
        left = right
        right *= 2
        number = base ** (- right)

    # left - x - right
    mid = None
    while left < right:
        mid = (left + right) // 2

        if base ** (-mid) + 1 != 1:
            left = mid + 1
        else:
            right = mid - 1

    return base ** (- right)


def assoc_test_additive():
    a = 1.0
    u = machine_precision_linear()

    b = c = u

    x1 = (a + b) + c
    x2 = a + (b + c)

    return x1, x2, x1 == x2


def assoc_test_multiplication(nr_tries=10000):
    for i in range(nr_tries):
        a = np.random.rand(3)
        if (a[0] * a[1]) * a[2] != a[0] * (a[1] * a[2]):
            return i, a, False

    return -1, -1, True


def test_multiplication_heuristic(nr_times):
    nr_tries = 10000

    average_index = 0
    min_index = nr_tries + 1
    max_index = -1
    found = 0
    for i in range(nr_times):
        index, vars, cond = assoc_test_multiplication(nr_tries=nr_tries)
        average_index += index
        min_index = min(min_index, index)
        max_index = max(max_index, index)
        if not cond:
            found += 1

    print(
        f"Statistics Multiplication Associativity:\n"
        f"Nr Tests = {nr_times}\n"
        f"Avg Index = {average_index / nr_times}\n"
        f"Min Index = {min_index}\n"
        f"Max Index = {max_index}\n"
        f"Accuracy = {found / nr_times * 100}\n"
    )


def tan_lentz_approximation(x, eps=1e-12):
    b_0 = 0
    b_1 = 1

    a_1 = x
    a_2 = -x * x

    f_0 = b_0
    if f_0 == 0:
        f_0 = eps

    C_0 = f_0
    D_0 = 0
    j = 1
    while True:
        D_1 = b_1 + a_1 * D_0
        if D_1 == 0:
            D_1 = eps

        C_1 = b_1 + a_1 / C_0
        if C_1 == 0:
            C_1 = eps

        D_1 = 1 / D_1
        Delta_1 = C_1 * D_1
        f_1 = Delta_1 * f_0

        if np.abs(Delta_1 - 1) < eps:
            break

        f_0 = f_1
        C_0 = C_1
        D_0 = D_1
        b_1 = b_1 + 2
        if j == 1:
            a_1 = a_2

        j += 1

    return f_1


def tan_polynomial_approximation(x):
    # Reduction from (-oo, oo) -> (-PI/2, PI/2)
    x += PI_2
    x %= PI
    x -= PI_2

    # From (-PI/2, PI/2) -> (0, PI/2)
    sign = 1
    if x < 0:
        sign = -1
        x = -x

    # From (0, PI/2) -> (0, PI/4)
    interval_change = False
    if PI_4 <= x < PI_2:
        x = PI_2 - x
        interval_change = True

    x_2 = x * x
    x_3 = x * x_2
    x_4 = x_2 * x_2
    x_6 = x_4 * x_2

    c = np.array([0.33333333333333333, 0.133333333333333333, 0.053968253968254, 0.0218694885361552])
    p = np.array([1, x_2, x_4, x_6])

    p_x = np.dot(c, p)

    result = x + p_x * x_3
    result *= sign

    if interval_change:
        return 1 / result
    else:
        return result


def apply_func_timer(x, f):
    start_time = time.time()

    f_x = []
    for val in x:
        f_x.append(f(val))
    f_x = np.array(f_x)

    stop_time = time.time()
    return f_x, stop_time - start_time


def tan_statistic(nr_tests=1e4):
    nr_tests = int(nr_tests)
    x = np.random.uniform(low=np.nextafter(-PI_2, PI_2), high=PI_2, size=nr_tests)
    real_tan = np.tan(x)

    poly_approx = apply_func_timer(x=x, f=tan_polynomial_approximation)
    print("Polynomial Timer: ", poly_approx[1])
    print("Polynomial Error Sum Measure: ", np.sum(abs(real_tan-poly_approx[0])))
    print("Polynomial Error Avg Measure: ", np.average(abs(real_tan-poly_approx[0])))

    print()

    lentz_approx = apply_func_timer(x=x, f=tan_lentz_approximation)
    print("Lentz Timer: ", lentz_approx[1])
    print("Lentz Error Sum Measure: ", np.sum(abs(real_tan - lentz_approx[0])))
    print("Lentz Error Avg Measure: ", np.average(abs(real_tan - lentz_approx[0])))


def homework_1_presentation():
    # T1 E1
    print("T1, E1:\n Machine Precision Linear u=", machine_precision_linear())
    print(" Machine Precision Binary u=", machine_precision_binary(), "\n")

    # T1 E2:
    x1, x2, cond = assoc_test_additive()
    print(f"T1, E2:\n X1 = {x1}, X2 = {x2}\n Assoc = {cond} \n")

    i, vars, cond = assoc_test_multiplication()
    print(f"T1, E2 Bonus:\n Vars = {vars}\n Assoc = {cond} \n")

    test_multiplication_heuristic(nr_times=10000)

    # T1 E3:
    print(f"T1, E3 Polynomial Approximation: ")
    # print(tan_polynomial_approximation(-25 * PI / 3))
    # print(tan_lentz_approximation(-25 * PI / 3))
    # print(np.tan(-25 * PI / 3))

    tan_statistic()

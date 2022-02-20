import numpy as np
from numpy.random import default_rng
from bisect import bisect_left
from sympy import *
import matplotlib.pyplot as plt


def generate_points(input_line, n=10000, scale_power=7):
    x0 = input_line[0]
    xn = input_line[1]
    f = input_line[2]

    scale = 10**scale_power

    rng = default_rng()
    numbers = rng.choice(scale-1, size=n-1, replace=False)  # [0, scale-1) UNIQUE NUMBERS
    numbers = numbers + 1                                   # (0, scale)
    numbers = numbers / scale                               # (0, 1)
    numbers = numbers * (xn-x0)                             # (0, x[n]-x[0])
    numbers = numbers + x0                                  # (x[0], x[n])

    xs = [x0, xn]
    xs.extend(numbers)
    xs.sort()

    ys = [f(x) for x in xs]
    return xs, ys


def polynomial_interpolation(xs, ys, xf, m=5):
    pass


def binary_search(a, x):
    i = bisect_left(a, x)
    if i != len(a):
        return i
    else:
        return -1


def calculate_derivative(xf, g):
    x = Symbol('x')
    y = g(x)
    yprime = y.diff(x)
    f = lambdify(x, yprime, 'numpy')
    return f(xf)


def calculate_derivative_2(x_f, g, scale_power=10):
    scale = 10 ** scale_power
    x = np.asarray([x_f, x_f+1/scale])
    dx = x[1] - x[0]
    y = np.array(list(map(g, x)))
    dydx = np.gradient(y, dx)
    return dydx[0]


def horner_scheme(a, x_f):
    n = len(a)
    d = a[0]

    for i in range(1, n):
        d = a[i] + d*x_f

    return d


def spline_interpolation(xs, ys, x_f, d_a):
    index = binary_search(xs, x_f)
    index -= 1

    h_i = xs[index+1]-xs[index]
    A = d_a

    A_n = -A + 2*(ys[1]-ys[0]) / (xs[1] -xs[0])
    for i in range(index):
        A_n = -A + 2*(ys[i+1]-ys[i]) / (xs[i+1] -xs[i])

        if i == index-1:
            break

        A = A_n

    Sf = (A_n - A) / h_i * (x_f - xs[index])**2 + A_n*(x_f-xs[index]) + ys[index]

    return Sf


def homework_6_presentation():
    inputs = [
        [1, 5, lambda x: x*x - 12*x + 30],
        [0, 1.5, lambda x: np.sin(x) - np.cos(x)],
        [0, 2, lambda x: 2*x**3 - 3*x + 15],
        [2, 100, lambda x: np.sin(x) * np.log(x) - (2 + np.cos(x)) ** np.sqrt(x)],
    ]

    i = 2

    xs, ys = generate_points(inputs[i], n=int(6e2))

    x_f = 1.2
    Sf = spline_interpolation(xs, ys, x_f=x_f, d_a=calculate_derivative_2(inputs[i][0], inputs[i][2]))

    print(f"Real value of f({x_f}):\t\t\t\t\t {inputs[i][2](x_f)}")
    print(f"Spline Interpolation value of f({x_f}):\t {Sf}")

    print("Spline Interpolation error:\t\t\t\t", abs(Sf - inputs[i][2](x_f)))

    Sfs, x_Sfs = [], []
    step = 0.1
    j = inputs[i][0] + step
    while j <= inputs[i][1]:
        Sf = spline_interpolation(xs, ys, x_f=j, d_a=calculate_derivative_2(inputs[i][0], inputs[i][2]))
        Sfs.append(Sf)
        x_Sfs.append(j)
        j += step

    plt.plot(xs, ys)
    plt.plot(x_Sfs, Sfs)
    plt.show()

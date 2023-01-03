#!/usr/bin/env python3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

time0 = 7 * 60 + 35

times = np.array([7 * 60 + 35, 8 * 60, 11 * 60 + 35, 15 * 60 + 44, 16 * 60 + 36, 17 * 60 + 42]) - time0
temps = np.array([54.5, 53.6, 45.8, 40.4, 39.4, 38.3])

C0 = 25

def f(t, a, b):
    return a * np.exp(-b * t) + C0

def linear_fit(xs, ys, c):
    ys = np.log(ys - c)
    K, A = np.polyfit(xs, ys, 1)
    return np.exp(A), -K

lin = linear_fit(times, temps, C0)

def go_linear():
    print(lin, "Time constant (hours):", 1 / lin[1] / 60)
    A, K = lin
    plt.scatter(times, temps)
    xs = np.linspace(0, max(times), 100)
    plt.plot(xs, f(xs, *lin))
    plt.show()
    
fit, fitcov = curve_fit(f, times, temps, maxfev = 1000)

def go():
    print(fit)
    plt.scatter(times, temps)
    xs = np.linspace(0, max(times), 100)
    plt.plot(xs, f(xs, *fit))
    plt.show()


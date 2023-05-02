import numpy as np
import matplotlib.pyplot as plt


def x_tramo_1(t, x0):
    return x0 + 70 * t


def x_tramo_2(t, x0):
    return (x0 - 36) * np.exp(5 * t) + 36


def x_tramo_3(t, x0):
    return x0 + 120 * t


def x_tramo(t, x0, toggle):
    x_tramos = [x_tramo_1, x_tramo_2, x_tramo_3]
    return x_tramos[toggle](t, x0)


def x_loop(x0, delta_t):
    limits = [50, 60, np.inf]
    t0 = 0
    toggle = 1
    if x0 <= 50:
        toggle = 0
    elif x0 >= 60:
        toggle = 2
    output = []
    for t in np.arange(0, 1, delta_t):
        output.append(x_tramo(t - t0, x0, toggle))
        if output[-1] >= limits[toggle]:
            x0 = limits[toggle]
            toggle += 1
            t0 = t
    return output, np.arange(0, 1, delta_t)


def pintar_curvas(x0=0, xf=100, delta_x=2.5, t0=0, tf=1, delta_t=0.0001):
    for x_0 in np.arange(x0, xf, delta_x):
        x, y = x_loop(x_0, delta_t)
        plt.plot(x, y, "blue")
    plt.xlim([x0, xf])
    plt.ylim([t0, tf])
    plt.xlabel("Eje x (km)")
    plt.ylabel("Eje t (h)")
    plt.title("Curvas características de θ en el plano x-t")
    plt.show()


pintar_curvas()

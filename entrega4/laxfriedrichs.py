import numpy as np
import matplotlib.pyplot as plt


def crear_matriz(N, M):
    return np.zeros([M + 1, N + 1])


def u0(x):
    return 0.2 + 0.5 * np.exp(-((x - 50) / 30) ** 2)


def condis_iniciales_t(t):
    return 0.2


def f(u):
    return 120 * u * (1 - u)


def representar_mallado(mallado, title, xlabel='Eje x', ylabel='Eje t', N=1000, M=1000, plot_show=False, x_inicial=0,
                        x_final=100, t_inicial=0, t_final=1):
    # Funcion que genera un mapa de calor para una matriz dada. Permite seleccionar el titulo, el nombre de cada eje, los valores relativos al número de puntos y extremos (para poder marcar correctamente los ejes), y tiene un toggle plot_show por si queremos que se ejecute plt.show() al acabar, cosa que no querremos al buscar que se formen varios plots a la vez. La función está extraída con ligeras modificaciones de la entrega anterior.
    fig, ax = plt.subplots()
    im = ax.imshow(mallado, aspect='auto', cmap='hot')
    cbar = ax.figure.colorbar(im, ax=ax)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    y_ticks = np.linspace(t_inicial, t_final, 6)
    ax.set_yticks(np.linspace(0, M, 6))
    ax.set_yticklabels(np.round(y_ticks, decimals=1))
    ax.invert_yaxis()
    x_ticks = np.linspace(x_inicial, x_final, 6)
    ax.set_xticks(np.linspace(0, N, 6))
    ax.set_xticklabels(np.round(x_ticks, decimals=0))
    if plot_show:
        plt.show()


def calcular_lax(delta_x, x0=0, xf=100, t0=0, tf=2.5):
    delta_t = delta_x / 120
    N = int((xf - x0) / delta_x)
    M = int((tf - t0) / delta_t)
    parametro_lambda = delta_t / delta_x
    u = crear_matriz(N, M)
    for n in range(N + 1):
        u[0, n] = u0(x0 + delta_x * n)
    for m in range(M + 1):
        u[m, 0] = condis_iniciales_t(t0 + delta_t * m)
        u[m, N] = condis_iniciales_t(t0 + delta_t * m)
    for m in range(1, M + 1):
        for n in range(1, N):
            u[m, n] = (u[m - 1, n + 1] + u[n - 1, n - 1] - parametro_lambda * (
                        f(u[m - 1, n + 1]) - f(u[m - 1, n - 1]))) / 2
    representar_mallado(u, title='u(x,t) Lax-Friedrichs', N=N, M=M, plot_show=True, x_inicial=x0, x_final=xf,
                        t_inicial=t0, t_final=tf)
    return u


calcular_lax(0.2)

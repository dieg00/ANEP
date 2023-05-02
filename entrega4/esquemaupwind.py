import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera


def crear_matriz(N, M):
    return np.zeros([M + 1, N + 1])


def vs(x):
    if x <= 50:
        return 70
    elif x >= 60:
        return 120
    else:
        return 70 + 50 * (1 - (60 - x) / 10)


# def theta_value(matriz_theta,n,m,parametro_lambda,x0=0,delta_x=0.01,t0=0,delta_t=0.01):
#    matriz_theta[n,m] = matriz_theta[n-1,m]-vs(x0+n*delta_t)*parametro_lambda*(matriz_theta[n-1,m]-matriz_theta[n-1,m-1])
#    return matriz_theta

def condis_iniciales_x(x):
    return vs(x) * np.exp(-(np.abs(x - 10) ** 2) / 25)


def condis_iniciales_t(t):
    return 7


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


def calcular_theta(N, M, x0=0, xf=100, t0=0, tf=1, title="θ(x,t) Modelo Simplificado"):
    matriz_theta = crear_matriz(N, M)
    delta_x = (xf - x0) / N
    delta_t = (tf - t0) / M
    parametro_lambda = delta_t / delta_x
    for n in range(N + 1):
        matriz_theta[0, n] = condis_iniciales_x(x0 + delta_x * n)
    for m in range(M + 1):
        matriz_theta[m, 0] = condis_iniciales_t(t0 + delta_t * m)
    for n in range(1, N + 1):
        for m in range(1, M + 1):
            matriz_theta[m, n] = matriz_theta[m - 1, n] - vs(x0 + delta_x * n) * parametro_lambda * (
                        matriz_theta[m - 1, n] - matriz_theta[m - 1, n - 1])
    representar_mallado(matriz_theta, N=N, M=M, title=title, plot_show=True)
    return matriz_theta


def calcular_u(matriz_theta, x0=0, xf=100, title=" Densidad u(x,t) Modelo Simplificado"):
    M, N = matriz_theta.shape
    delta_x = (xf - x0) / N
    for n in range(N):
        matriz_theta[:, n] /= vs(x0 + n * delta_x)
    representar_mallado(matriz_theta, N=N, M=M, title=title, plot_show=True)
    return matriz_theta


def animar(matriz_datos, x0=0, xf=100):
    M, N = matriz_datos.shape
    x = np.linspace(x0, xf, N)
    fig = plt.figure()
    camera = Camera(fig)
    plt.xlim(x0, xf)
    for m in range(M):
        if m % 10:
            plt.plot(x, matriz_datos[m, :], "b", linewidth=5)
            plt.xlabel('Eje x(km)')
            plt.ylabel('Densidad de tráfico u')
            camera.snap()
    animation = camera.animate(interval=1, blit=False)
    animation.save('animacion1.gif', writer='Pillow')  # esto guarda el video
    return animation


theta = calcular_theta(1000, 1200)
calcular_theta(1000, 1000, title="θ(x,t) violando CFL")
u = calcular_u(theta)
animar(u)

import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera


def crear_matriz(N, M):
    # Crea un array vacío que servirá como array en el que guardar los valores de las funciones en los puntos del mallado.

    return np.zeros([M + 1, N + 1])


def u0(x):
    # Define la función u0 de las condiciones espaciales de frontera, que son las condiciones iniciales

    return 0.2 + 0.5 * np.exp(-((x - 50) / 30) ** 2)


def condis_frontera_t(t):
    # Al igual que en esquemaupwind.py, podriamos haber simplemente asignado la constante, pero se ha optado por dejar la función de cara a posibles modificaciones que se quiera hacer

    return 0.2


def f(u):
    # Define la función f que servirá para iterar en el método Lax-Friedrichs

    return 120 * u * (1 - u)


def representar_mallado(mallado, title, xlabel="Eje x (km)", ylabel="Eje t (h)", N=1000, M=1000, plot_show=False,
                        x_inicial=0,
                        x_final=100, t_inicial=0, t_final=1):
    # Funcion que genera un mapa de calor para una matriz dada. Permite seleccionar el titulo, el nombre de cada eje, los valores relativos al número de puntos y extremos (para poder marcar correctamente los ejes), y tiene un toggle plot_show por si queremos que se ejecute plt.show() al acabar, cosa que no querremos al buscar que se formen varios plots a la vez. La función está extraída con ligeras modificaciones de la entrega anterior.

    fig, ax = plt.subplots()
    im = ax.imshow(mallado, aspect="auto", cmap="hot")
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


def calcular_lax(delta_x, x0=0, xf=100, t0=0, tf=2.5, title="Densidad u(x,t) No Lineal, Lax-Friedrichs"):
    # Calcula la solución u del modelo no lineal a partir del método de Lax-Friedrichs, para ello requiere definir el intervalo del mallado espacial. También se pueden ajustar los límites espaciales y temporales de la región a estudio, así como el título de la gráfica generada.

    # Con delta_x podemos calcular delta_t por requerimientos del enunciado, así como las dimensiones del mallado. En consecuencia, generamos la matriz que contendrá los valores de dicho mallado.
    delta_t = delta_x / 120
    N = int((xf - x0) / delta_x)
    M = int((tf - t0) / delta_t)
    u = crear_matriz(N, M)

    # Establecemos las condiciones de frontera
    for n in range(N + 1):
        u[0, n] = u0(x0 + delta_x * n)
    for m in range(M + 1):
        u[m, 0] = condis_frontera_t(t0 + delta_t * m)

        # Cabe destacar que aquí la condición se establecerá en el extremo derecho al igual que el izquierdo, a diferencia de lo que sucedía en esquemaupwind.py
        u[m, N] = condis_frontera_t(t0 + delta_t * m)

    # Iteramos con el método de Lax-Friedrichs para ir populando la matriz con los valores de la solución
    for m in range(1, M + 1):
        for n in range(1, N):
            u[m, n] = (u[m - 1, n + 1] + u[m - 1, n - 1] - (delta_t / delta_x) * (f(u[m - 1, n + 1]) - f(u[m - 1, n - 1]))) / 2
    representar_mallado(u, title=title, N=N, M=M, plot_show=True, x_inicial=x0, x_final=xf, t_inicial=t0, t_final=tf)
    return u


def animar(matriz_datos, x0=0, xf=100, y0=0, yf=1, t0=0, tf=2.5):
    # Copiada de esquemaupwind.py, permite animar temporalmente la solución. Se puede encontrar en dicho archivo una explicación de su funcionamiento

    M, N = matriz_datos.shape
    delta_t = (tf - t0) / M
    x = np.linspace(x0, xf, N)
    fig = plt.figure()
    camera = Camera(fig)
    plt.xlim([x0, xf])
    plt.ylim([y0, yf])
    plt.xlabel("Eje x (km)")
    plt.ylabel("Densidad u")
    plt.title("Evolución temporal u(x,t) No Lineal")
    for m in range(M)[::25]:
        plt.plot(x, matriz_datos[m, :], "b", linewidth=5)
        plt.text((xf - x0) * 0.55 + x0, (yf - y0) * 0.9 + y0, f"Tiempo transcurrido: {m * delta_t:.2f} horas",
                 fontsize=8.5)
        camera.snap()
    animation = camera.animate(interval=1, blit=False)
    animation.save("animacion2.gif", writer="Pillow")
    return animation


# Calculamos la solución u para el delta_x dado en el enunciado, y posteriormente se anima
u = calcular_lax(0.2)
animar(u)

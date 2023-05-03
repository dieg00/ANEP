import numpy as np
import matplotlib.pyplot as plt
from celluloid import Camera


def crear_matriz(N, M):
    # Crea un array vacío que servirá como array en el que guardar los valores de las funciones en los puntos del mallado.

    return np.zeros([M + 1, N + 1])


def vs(x):
    # Define la velocidad de los coches en un punto determinado, de acuerdo con las condiciones de límite de velocidad de la vía.

    if x <= 50:
        return 70
    elif x >= 60:
        return 120
    else:
        return 70 + 50 * (1 - (60 - x) / 10)


def condis_frontera_x(x):
    # Define las condiciones de frontera espaciales, que son las condiciones iniciales.

    return vs(x) * np.exp(-(np.abs(x - 10) ** 2) / 25)


def condis_frontera_t(t):
    # Define las condiciones de frontera temporales, que se encuentran en el extremo izquierdo. Aunque en este caso sea un valor constante, se ha preferido definir la función por modularidad en caso de creer modificarlas a posteriori.

    return 7


def representar_mallado(mallado, title, xlabel='Eje x (km)', ylabel='Eje t (h)', N=1000, M=1200, plot_show=False,
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


def calcular_theta(N, M, x0=0, xf=100, t0=0, tf=1, title="θ(x,t) Modelo Simplificado"):
    # Función que permite calcular la matriz theta solución de su correspondiente problema. Requiere especificar el número de divisiones espaciales y temporales que se van a requerir, así como los extremos de la región y el título de la gráfica que se dibuja.

    # Creamos la matriz y calculamos los intervalos
    matriz_theta = crear_matriz(N, M)
    delta_x = (xf - x0) / N
    delta_t = (tf - t0) / M

    # Aplicamos las condiciones de frontera espaciales y temporales
    for n in range(N + 1):
        matriz_theta[0, n] = condis_frontera_x(x0 + delta_x * n)
    for m in range(M + 1):
        matriz_theta[m, 0] = condis_frontera_t(t0 + delta_t * m)

    # Iteramos siguiendo el esquema upwind comentado en la práctica.
    for n in range(1, N + 1):
        for m in range(1, M + 1):
            matriz_theta[m, n] = matriz_theta[m - 1, n] - vs(x0 + delta_x * n) * (delta_t / delta_x) * (
                    matriz_theta[m - 1, n] - matriz_theta[m - 1, n - 1])

    # Representamos la solución en un mapa de calor
    representar_mallado(matriz_theta, N=N, M=M, title=title, plot_show=True)

    # El output es la matriz con los valores, que se usará posteriormente
    return matriz_theta


def calcular_u(matriz_theta, x0=0, xf=100, title="Densidad u(x,t) Modelo Simplificado"):
    # Utilizando la matriz theta obtenida con calcular_theta, permite deshacer el cambio de variable y graficar la nueva matriz u

    # Necesitamos las dimensiones de la matriz y el tamaño del intervalo espacial, de cara a calcular vs(x)
    M, N = matriz_theta.shape
    delta_x = (xf - x0) / N

    # Dividimos entre vs(x) para deshacer el cambio de variable
    for n in range(N):
        matriz_theta[:, n] /= vs(x0 + n * delta_x)

    # Graficamos la nueva matriz, que a pesar de mantener el nombre de la variable, se trata ya de u
    representar_mallado(matriz_theta, N=N, M=M, title=title, plot_show=True)

    # Como output tenemos la matriz u, que se usará para animar la situación
    return matriz_theta


def animar(matriz_u, x0=0, xf=100, y0=0, yf=1, t0=0, tf=1):
    # Esta función permite animar la solución u en el tiempo, para poder observar mejor la evolución de dicha solución. Necesitamos la propia matriz con la solución u, así como los extremos espaciales (x0,xf), temporales (t0,tf) de la región a estudio, y los extremos (y0,yf) de la función u, de cara a la gráfica.

    # Necesitamos conocer las dimensiones del mallado
    M, N = matriz_u.shape

    # Obtenemos la división del intervalo espacial de acuerdo a dicho mallado
    x = np.linspace(x0, xf, N)

    # Calculamos el tamaño de división del intervalo temporal
    delta_t = (tf - t0) / M

    # Establecemos los valores constantes del plot que se animará
    fig = plt.figure()
    camera = Camera(fig)
    plt.xlim([x0, xf])
    plt.ylim([y0, yf])
    plt.xlabel("Eje x (km)")
    plt.ylabel("Densidad u")
    plt.title("Evolución temporal u(x,t) Simplificado")

    # Iteramos en slices temporales para ir viendo la evolución, nótese que para no generar un gif demasiado largo sólo tomamos un "frame" de cada 25.
    for m in range(M)[::25]:
        # Graficamos el slice temporal
        plt.plot(x, matriz_u[m, :], "b", linewidth=5)

        # Iremos modificando un cuadro de texto que permita ver el tiempo transcurrido en cada frame, para poder localizar temporalmente los momentos interesantes
        plt.text((xf - x0) * 0.55 + x0, (yf - y0) * 0.9 + y0, f"Tiempo transcurrido: {m * delta_t:.2f} horas",
                 fontsize=8.5)

        # Captura un frame
        camera.snap()

    # Animamos y guardamos la animación
    animation = camera.animate(interval=1, blit=False)
    animation.save("animacion1.gif", writer="Pillow")


# Calculamos la solución theta si se cumple CFL
theta = calcular_theta(1000, 1200)

# Probamos a violar CFL
calcular_theta(1000, 1000, title="θ(x,t) violando CFL")

# Con la theta correctamente calculada obtenemos u, y posteriormente la animamos
u = calcular_u(theta)
animar(u)

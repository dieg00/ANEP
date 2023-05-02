import numpy as np
import matplotlib.pyplot as plt


def x_tramo_1(t, x0):
    # Funcion x(t) en el tramo [0,50]
    return x0 + 70 * t


def x_tramo_2(t, x0):
    # Función x(t) en el tramo [50,60]
    return (x0 - 36) * np.exp(5 * t) + 36


def x_tramo_3(t, x0):
    # Función x(t) en el tramo [60,100]
    return x0 + 120 * t


def x_tramo(t, x0, toggle):
    # Con un toggle [0,1,2] se elige el tramo actual, permitiendo unir las tres funciones de tramos en una
    x_tramos = [x_tramo_1, x_tramo_2, x_tramo_3]
    return x_tramos[toggle](t, x0)


def calcular_curva(x0, delta_t):
    # Calcula las coordenadas (x,t) de la curva característica que tiene a x0 como valor inicial.

    #Definimos los límites de los intervalos
    limits = [50, 60, np.inf]

    # El tiempo inicial en una x dada será el instante en que se haya alcanzado el valor inicial de su intervalo o 0 si se partía originalmente de dicho intervalo.
    t0 = 0
    toggle = 1
    if x0 <= 50:
        toggle = 0
    elif x0 >= 60:
        toggle = 2

    # Aquí se guardarán los valores a calcular
    output_x = []

    # Para cada instante asignamos el valor correspondiente según la función en dicho tramo. Si superamos el límite del tramo, entonces x0 pasa a ser el valor inicial del nuevo tramo, el toggle de la función cambia y t0 pasa a ser dicho instante en que se ha cambiado de tramo.
    for t in np.arange(0, 1, delta_t):
        output_x.append(x_tramo(t - t0, x0, toggle))
        if output_x[-1] >= limits[toggle]:
            x0 = limits[toggle]
            toggle += 1
            t0 = t

    # Devolvemos tanto el vector con los valores de x como el de t, de cara a graficar la curva
    return output_x, np.arange(0, 1, delta_t)


def pintar_curvas(x0=0, xf=100, delta_x=2.5, t0=0, tf=1, delta_t=0.0001):
    # Esta función va iterando por los distintos valores de posiciones iniciales x_0 que puede haber y dibujando la curva característica de cada uno.
    for x_0 in np.arange(x0, xf, delta_x):
        x, y = calcular_curva(x_0, delta_t)
        plt.plot(x, y, "blue")
    plt.xlim([x0, xf])
    plt.ylim([t0, tf])
    plt.xlabel("Eje x (km)")
    plt.ylabel("Eje t (h)")
    plt.title("Curvas características de θ en el plano x-t")
    plt.show()


pintar_curvas()

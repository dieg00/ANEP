import numpy as np
import matplotlib.pyplot as plt


def u0(x0):
    # Condiciones de frontera espaciales (iniciales)
    return 0.2 + 0.5 * np.exp(-((x0 - 50) / 30) ** 2)


def x_curva(t, x0, u0_valor):
    # Me calcula la posicion a lo largo de una curva caracteristica, de acuerdo a la expresión x=x0+f'(u0)·t
    return x0 + 120 * (1 - 2 * u0_valor) * t


def calcular_curva(x0, t0=0, tf=1, M=1000):
    # Itera para una posición inicial dada valores de tiempo para ir obteniendo un vector de posiciones con el que graficar la curva, de forma a análoga a como funcionaba en curvascaracteristicasB.py

    u0_valor = u0(x0)
    output_x = np.zeros(M)
    tiempos = np.linspace(t0, tf, M)
    for m in range(M):
        output_x[m] = x_curva(tiempos[m], x0=x0, u0_valor=u0_valor)
    return output_x, tiempos


def pintar_curvas(t0=0, tf=1, M=1000, x0=0, xf=100, delta_x=2.5):
    # Grafica las curvas caracteristicas iterando entre distintos valores de la condición inicial

    for x_0 in np.arange(x0, xf, delta_x):
        x, y = calcular_curva(x_0, t0=t0, tf=tf, M=M)
        plt.plot(x, y, "blue")
    plt.xlim([x0, xf])
    plt.ylim([t0, tf])
    plt.xlabel("Eje x (km)")
    plt.ylabel("Eje t (h)")
    plt.title("Curvas características de u en el plano x-t")
    plt.show()


pintar_curvas()

import numpy as np
import matplotlib.pyplot as plt


def f(x):  # Función de la que calcular su derivada
    return np.cos(x)


def Df(x):  # Derivada exacta, obtenida de forma teórica
    return -np.sin(x)


def Delta_forward_h(f, x, h):  # Diferencias finitas de forma forward
    return (-11 / 6 * f(x) + 3 * f(x + h) - 3 / 2 * f(x + 2 * h) + 1 / 3 * f(x + 3 * h)) / h


def delta_centrada_h(f, x, h):  # Diferencias finitas centradas
    return (1 / 12 * f(x - 2 * h) - 2 / 3 * f(x - h) + 2 / 3 * f(x + h) - 1 / 12 * f(x + 2 * h)) / h


def error(f, derivada, aprox, h, intervalo=[0, 2 * np.pi], pasos=2000):
    # Toma la función f, su derivada exacta, el metodo 'aprox' de diferencias finitas, la h a usar,
    # así como el intervalo de estudio y el número de pasos, y calcula el error en la aproximación
    x = np.linspace(intervalo[0], intervalo[1], pasos) #Sampleamos unos x equiespaciados en el intervalo
    return np.sum(abs(derivada(x) - aprox(f, x, h))) #Obtenemos el tamaño del vector de error (valor absoluto de la
    # diferencia entre la derivada y su aproximación por diferencias finitas) por medio de la métrica infinito.


h_array = np.logspace(-3, -1, 1000) #Equiespaciamos logaritmicamente los valores de h
erroresD = np.zeros(len(h_array))
erroresd = np.zeros(len(h_array))

for i in range(len(h_array)):
    erroresD[i] = error(f, Df, Delta_forward_h, h_array[i])#Usando dif forward
    erroresd[i] = error(f, Df, delta_centrada_h, h_array[i]) #Usando dif centrada

#Graficamos con ejes logarítmicos para apreciar la diferencia entre las pendientes
plt.loglog(h_array, erroresD, 'r', label='Diferencias forward')
plt.loglog(h_array, erroresd, 'b', label='Diferencias centradas')
plt.xlabel('Diferencia finita h')
plt.ylabel('Error de aproximación')
plt.legend()
plt.title('Errores según la diferencia finita')
plt.show()

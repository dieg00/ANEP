import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.cos(x)


def derivada_exacta(x):
    return -np.sin(x)


def Delta_h(f,x,h):
    return (-11/6*f(x)+3*f(x+h)-3/2*f(x+2*h)+1/3*f(x+3*h))/h


def delta_h(f,x,h):
    return (1/12*f(x-2*h)-2/3*f(x-h)+2/3*f(x+h)-1/12*f(x+2*h))/h


def error(f,derivada,aprox,h,intervalo=[0, 2*np.pi],pasos=400):
    x=np.linspace(intervalo[0],intervalo[1], pasos)
    return np.sum(abs(derivada(x)-aprox(f,x,h)))

print(error(f,derivada_exacta,Delta_h,1))
h=np.linspace(0.001,1,1000)
#plt.plot(h,error(f,Delta_h,h))

plt.show()


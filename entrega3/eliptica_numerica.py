import numpy as np


def f(x,y,C=40):
    if np.sqrt(x**2+y**2)<=1/2:
        return C*(1-2*np.sqrt(x**2+y**2))
    else:
        return 0

def a_i(xi,k):
    if np.abs(xi) <= 1:
        return 1
    else:
        return np.exp(-k*(np.abs(xi)-1)**2)

def a(x,y,k_x=25):
    return a_i(x,k_x)

def b(x,y,k_y=25):
    return a_i(y,k_y)

def generar_mallado(funcion,x_inicial,x_final,y_inicial,y_final,N,M):
    output = np.zeros()
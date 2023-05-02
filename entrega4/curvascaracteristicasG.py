import numpy as np
import matplotlib.pyplot as plt

def u0(x0):
    return 0.2+0.5*np.exp(-((x0-50)/30)**2)

def x(t,x0,u0_valor):
    return x0+120*(1-2*u0_valor)*t

def calcular_curva(x0,t0=0,tf=1,M=1000):
    u0_valor=u0(x0)
    output = np.zeros(M)
    tiempos = np.linspace(t0,tf,M)
    for m in range(M):
        output[m] = x(tiempos[m],x0=x0,u0_valor=u0_valor)
    return output, tiempos

def main_loop(t0=0,tf=1,M=1000, x0=0,xf=100, delta_x=2.5):
    for x0 in np.arange(x0,xf,delta_x):
        x,y = calcular_curva(x0,t0=t0,tf=tf,M=M)
        plt.plot(x, y, "blue")
    plt.xlim([0, 100])
    plt.ylim([0, 1])
    plt.xlabel("Dimensión espacial x (km)")
    plt.ylabel("Dimensión temporal t (h)")
    plt.title("Curvas características en el plano x-t")
    plt.show()


main_loop()
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import random as rand


def f(x,y,C=40):
    if np.sqrt(x**2+y**2)<=1/2:
        return C*(1-2*np.sqrt(x**2+y**2))
    else:
        return 0

def a_i(xi,k=25):
    if np.abs(xi) <= 1:
        return 1
    else:
        return np.exp(-k*(np.abs(xi)-1)**2)

def a(x,y,k_x=25):
    return a_i(x,k_x)

def b(x,y,k_y=25):
    return a_i(y,k_y)

def d(x,y):
    return x
def generar_mallado(funcion,x_inicial=-3/2,x_final=3/2,y_inicial=-3/2,y_final=3/2,N=1000,M=2000):
    output = np.zeros([M, N])
    h1 = (x_final-x_inicial)/(N+1)
    h2 = (y_final-y_inicial)/(M+1)
    for i in range(N):
        for j in range(M):
            output[j,i]=funcion(x_inicial+i*h1,y_inicial+j*h2)
    return output

def representar_mallado(mallado,title,xlabel='Eje X',ylabel='Eje Y',N=1000,M=2000):
    # Funcion que genera el mapa de calor, para un sistema, y denotando unas condiciones
    fig, ax = plt.subplots()
    im = ax.imshow(mallado, aspect='auto', cmap='hot')
    cbar = ax.figure.colorbar(im, ax=ax)
    #ax.set_title("Sistema u(x,t), condiciones ")
    #ax.set_ylabel("Eje espacial (x)")
    #ax.set_xlabel("Eje temporal (t)")
    ax.set_title(title)
    ax.set_ylabel(xlabel)
    ax.set_xlabel(ylabel)
    y_ticks = np.linspace(-3/2, 3/2, 11)
    ax.set_yticks(np.linspace(0, M, 11))
    ax.set_yticklabels(np.round(y_ticks, decimals=1))
    x_ticks = np.linspace(-3/2, 3/2, 11)
    ax.set_xticks(np.linspace(0, N, 11))
    ax.set_xticklabels(np.round(x_ticks, decimals=1))
    #plt.show()


def T(func,N,xi_inicial=-3/2,xi_final=3/2,k=25):
    hi = (xi_final-xi_inicial)/(N+1)
    x_vector = np.linspace(xi_inicial+hi/2,xi_final-hi/2,N+1)
    vector_func = [0]*(N+1)
    for i in range(len(x_vector)):
        vector_func[i] = func(x_vector[i],k=25)
    output = (1/hi**2)*(np.diag(vector_func[:-1])+np.diag(vector_func[1:])-np.diag(vector_func[1:-1],-1)-np.diag(vector_func[1:-1],1))
    return output

def calcular_U(C=40,kx=25,ky=25,N=1000,M=2000,x_inicial=-3/2,x_final=3/2,y_inicial=-3/2,y_final=3/2):
    F = generar_mallado(f,x_inicial,x_final,y_inicial,y_final,N,M)
    Tx = T(a_i,N,x_inicial,x_final,kx)
    Ty = T(a_i,M,y_inicial,y_final,ky)
    U = sp.linalg.solve_sylvester(Ty,Tx,F)
    return U

f_fig={
    "func": f,
    "title": "f(x,y)"
}
a_data={
    "func": a,
    "title": "a(x,y)"
}
b_data={
    "func": b,
    "title": "b(x,y)"
}

def representar_solucion(C=40,kx=25,ky=25,N=1000,M=2000,x_inicial=-3/2,x_final=3/2,y_inicial=-3/2,y_final=3/2, titulo='Solución U(x,y)'):
    U = calcular_U(C,kx,ky,N,M,x_inicial,x_final,y_inicial,y_final)
    representar_mallado(U,titulo,N=N,M=M)


representar_solucion()


plt.show()
#print(T(a_i,5))
#for funcion in [f_fig,a_fig,b_fig]:
    #representar_mallado(generar_mallado(funcion["func"]),title=funcion["title"])
#plt.show()


#print(generar_mallado(d))
#print(generar_mallado(f))

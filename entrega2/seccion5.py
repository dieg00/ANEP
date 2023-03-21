import numpy as np
import matplotlib.pyplot as plt


def calcular_matriz(parametro_lambda, parametro_theta,N,Neumann=True,devolver_inversa=True):
    multiplicador = 1+Neumann
    N = ajuste_N(N,Neumann)
    matriz = np.diag([1+2*parametro_lambda*parametro_theta]*N)+np.diag([-parametro_theta*parametro_lambda]*(N-1),-1)+np.diag([-parametro_theta*parametro_lambda]*(N-1),1)
    matriz[0,1] *= multiplicador
    matriz[N-1,N-2] *= multiplicador
    if devolver_inversa:
        matriz = np.linalg.inv(matriz)
    return matriz

def generar_mallado(N,M,Neumann=True):
    mallado = np.zeros([N+1,M+1])
    return mallado

def generar_valores(funcion,x_inicial,x_final,N,Neumann=True):
    xi = np.linspace(x_inicial,x_final,N+1)
    return funcion(xi)

def inicializar_mallado(condiciones_iniciales,N,M,Neumann=True):
    u = generar_mallado(N=N,M=M)
    u[:,0] = condiciones_iniciales
    return u


def calcular_bi(u,i,j,f,N,parametro_lambda,parametro_theta, Delta_t, Neumann = True):
    if Neumann and i in [0,N]:
        i0 = min(i + 1, N - 1)
        termino_variable = 2 * u[i0, j]
    else:
        termino_variable = u[i+1,j]+u[i-1,j]
    return parametro_lambda*(1-parametro_theta)*termino_variable+u[i,j]*(1-2*parametro_lambda*(1-parametro_theta))+Delta_t*(parametro_theta*f[i,j+1]+parametro_theta*f[i,j])

def ajuste_N(N,Neumann=True):
    if Neumann:
        return N+1
    else:
        return N-1
def calcular_b(u,j,f,N,parametro_lambda,parametro_theta, Delta_t, Neumann = True):
    b = np.zeros([ajuste_N(N,Neumann),1])
    for i in range(N+1):
        if Neumann or i not in [0,N]:
            b[i-(not Neumann)]=calcular_bi(u=u,i=i,j=j,f=f,N=N,parametro_lambda=parametro_lambda,parametro_theta=parametro_theta,Neumann=Neumann, Delta_t=Delta_t)
    return b

def iterar_operaciones(u,matriz,N,M,f,parametro_lambda,parametro_theta,Delta_t,Neumann=True):
    #N = ajuste_N(N,Neumann)
    for j in range(M):
        bj = calcular_b(u=u,j=j,f=f,N=N,parametro_lambda=parametro_lambda,parametro_theta=parametro_theta,Delta_t=Delta_t,Neumann=Neumann)
        u[not Neumann:N+Neumann,j+1]=matriz.dot(bj).ravel()
    return u

def total_wrap(parametro_theta,N,M,T,x0=0,xf=1,T0=0,Neumann=True,condiciones_son_funcion=True,funcion_condiciones=np.sin, condiciones_iniciales=np.nan):
    Delta_t = (T-T0)/M
    Delta_x = (xf-x0)/N
    parametro_lambda = Delta_t/(Delta_x**2)
    if (1-parametro_theta)*parametro_lambda <= 1/2:
        raise ValueError("El parametro lambda no cumple las condiciones de convergencia")
    if condiciones_son_funcion:
        condiciones_iniciales=generar_valores(funcion=funcion_condiciones,x_inicial=x0,x_final=xf,N=N,Neumann=Neumann)
    u = inicializar_mallado(condiciones_iniciales=condiciones_iniciales,N=N,M=M,Neumann=Neumann)
    matriz = calcular_matriz(parametro_lambda=parametro_lambda, parametro_theta=parametro_theta,N=N,Neumann=Neumann,devolver_inversa=True)
    u = iterar_operaciones(u, matriz, N, M, np.zeros_like(u), parametro_lambda, parametro_theta, Delta_t, Neumann=Neumann)
    print(u.shape)
    return u

total_wrap(parametro_theta=0,N=10,M=5,T=5,funcion_condiciones=np.cos,Neumann=False)


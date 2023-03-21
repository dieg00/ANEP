import numpy as np
import matplotlib.pyplot as plt


def matriz(parametro_lambda, parametro_theta,N,Neumann=True,devolver_inversa=True):
    if Neumann:
        N+=1
        multiplicador=2
    else:
        multiplicador=1
        N-=1
    matriz = np.diag([1+2*parametro_lambda*parametro_theta]*N)+np.diag([-parametro_theta*parametro_lambda]*(N-1),-1)+np.diag([-parametro_theta*parametro_lambda]*(N-1),1)
    matriz[0,1] *= multiplicador
    matriz[N-1,N-2] *= multiplicador
    if devolver_inversa:
        matriz = np.linalg.inv(matriz)
    return matriz

def generar_mallado(N,M,Neumann=True):
    mallado = np.zeros([N+1,M+1])
    return mallado

def generar_condiciones_iniciales(funcion_inicial,x_inicial,x_final,N):
    xi = np.linspace(x_inicial,x_final,N+1)
    return funcion_inicial(xi)

def inicializar_mallado(condiciones_iniciales,N,M,Neumann=True):
    u = generar_mallado(N,M,Neumann)
    u[:,0] = condiciones_iniciales
    return u

def calcular_b():
    return


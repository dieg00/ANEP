#Importamos las librerias necesarias
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sp
import random as rand


def f(x,y,C=40):
    # Funcion f(x,y) con su parametro C
    if np.sqrt(x**2+y**2)<=1/2:
        return C*(1-2*np.sqrt(x**2+y**2))
    else:
        return 0


def a_i(xi,k=25):
    # Aprovechamos que tanto a como b son idénticas salvo la variable y la constante k, para crear una funcion que calcule la formula con una sola variable y la k.
    if np.abs(xi) <= 1:
        return 1
    else:
        return np.exp(-k*(np.abs(xi)-1)**2)


# Entendiendo a y b como funciones de dos variables (para poder crear la funcion representar_mallado que sirva para ambas y llame a cada una dándole valores de x e y), que en el fondo solo llaman a la funcion a_i con la variable y k correspondientes.
def a(x,y,k_x=25): 
    # Funcion a(x,y)=a(x) con su constante k_x
    return a_i(x,k_x)

def b(x,y,k_y=25):
    # Funcion b(x,y)=b(y) con su constante k_y
    return a_i(y,k_y)


def generar_mallado(funcion,parametro,x_inicial=-3/2,x_final=3/2,y_inicial=-3/2,y_final=3/2,N=1000,M=2000):
    # La funcion generar mallado crea una matriz con los valores de una cierta funcion en una malla bidimensional. Basta darle la funcion, las coordenadas de la frontera en cada eje, y el numero de puntos interiores que queremos tenga la malla (los de la frontera para nuestro problema serán nulos y no será necesario tratarlos)

    # Generamos la matriz
    output = np.zeros([M, N])

    # Calculamos la distancia en cada eje entre puntos consecutivos de la malla. Nótese que al incluir los dos puntos de frontera para una línea pasamos a tener N+2 puntos con N+1 intervalos
    h1 = (x_final-x_inicial)/(N+1)
    h2 = (y_final-y_inicial)/(M+1)

    # Poblamos la matriz con el valor en cada punto de la malla
    for i in range(N):
        for j in range(M):
            output[j,i]=funcion(x_inicial+i*h1,y_inicial+j*h2, parametro)
    return output


def representar_mallado(mallado,title,xlabel='Eje X',ylabel='Eje Y',N=1000,M=2000, plot_show=False, x_inicial=-3/2, x_final=3/2, y_inicial=-3/2, y_final=3/2):
    # Funcion que genera un mapa de calor para una matriz dada. Permite seleccionar el titulo, el nombre de cada eje, los valores relativos al número de puntos y extremos (para poder marcar correctamente los ejes), y tiene un toggle plot_show por si queremos que se ejecute plt.show() al acabar, cosa que no querremos al buscar que se formen varios plots a la vez. La función está extraída con ligeras modificaciones de la entrega anterior.
    fig, ax = plt.subplots()
    im = ax.imshow(mallado, aspect='auto', cmap='hot')
    cbar = ax.figure.colorbar(im, ax=ax)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    y_ticks = np.linspace(y_inicial, y_final, 11)
    ax.set_yticks(np.linspace(0, M, 11))
    ax.set_yticklabels(np.round(y_ticks, decimals=1))
    x_ticks = np.linspace(x_inicial, x_final, 11)
    ax.set_xticks(np.linspace(0, N, 11))
    ax.set_xticklabels(np.round(x_ticks, decimals=1))
    if plot_show:
        plt.show()


def T(N,xi_inicial=-3/2,xi_final=3/2,func=a_i,k=25):
    # Funcion que permite calcular ambas matrices Tx y Ty para la resolución por el método de sylvester. Requiere como input el numero de puntos interiores en el eje correspondiente, la funcion que se utiliza (al ser las a o b se puede introducir la a_i seleccionando la coordenada y k correspondientes), los limites en el eje y la constante k.

    # Computamos el intervalo entre puntos de la malla del eje correspondiente
    hi = (xi_final-xi_inicial)/(N+1)

    # Valores de las coordenadas en los que evaluar a_i para obtener el vector correspondiente a_k o b_k. Ver abajo explicación de los límites
    x_vector = np.linspace(xi_inicial+hi/2,xi_final-hi/2,N+1)
    
    # Creamos el vector que contenga los valores discretos de la funcion
    vector_func = [0]*(N+1)

    # Poblamos el vector con los valores, que van desde a_i_{1/2} hasta a_i_{N+1/2} con intervalo hi, es decir, hi/2 antes y después de las coordenadas en las frontera
    for i in range(len(x_vector)):
        vector_func[i] = func(xi=x_vector[i],k=k)

    # Aplicamos la fórmula hallada en teoría, aprovechando que salvo seleccionar el vector de la funcion a_i apropiada (y el paso hi), la fórmula es idéntica
    output = (1/hi**2)*(np.diag(vector_func[:-1])+np.diag(vector_func[1:])-np.diag(vector_func[1:-1],-1)-np.diag(vector_func[1:-1],1))

    # Devuelve la matriz
    return output


def calcular_U(C=40,kx=25,ky=25,N=1000,M=2000,x_inicial=-3/2,x_final=3/2,y_inicial=-3/2,y_final=3/2):
    # Funcion que permite calcular la solucion U atendiendo al mallado seleccionado y las constantes involucradas

    # Generamos la matriz F que no es nada mas que la evaluacion de f(x,y) en el mallado
    F = generar_mallado(f,C,x_inicial,x_final,y_inicial,y_final,N,M)

    # Calculamos las matrices Tx y Ty utilizando la funcion T que definimos ambas que aprovechaba la similitud entre ambas
    Tx = T(N,x_inicial,x_final,k=kx)
    Ty = T(M,y_inicial,y_final,k=ky)

    # Usamos el método de sylvester programado en la libreria scipy para calcular la matriz U
    U = sp.solve_sylvester(Ty,Tx,F)

    #Recuperamos dicha matriz U
    return U


def representar_solucion(C=40,kx=25,ky=25,N=1000,M=2000,x_inicial=-3/2,x_final=3/2,y_inicial=-3/2,y_final=3/2, titulo='Solución U(x,y)', plot_show=False):
    # Esta funcion permite de forma directa computar y representar la solución U dado un mallado y condiciones.

    # Calculamos primero U
    U = calcular_U(C,kx,ky,N,M,x_inicial,x_final,y_inicial,y_final)

    # Con el U calculado lo representamos graficamente. Hemos añadido el titulo como posible variable para cuando queramos iterar varias posibles consonantes.
    representar_mallado(U,titulo,N=N,M=M,plot_show=plot_show)

    # La funcion devuelve la matriz solucion
    return U


# Toggles para elegir qué apartados se ejecutan
apartado_h=0
apartado_i=0
apartado_j=0
apartado_k=1


if apartado_h:

    # Utilizamos los diccionarios para simplificar el for loop
    f_data={
        "func": f,
        "title": "f(x,y)",
        "par": 40,
    }
    a_data={
        "func": a,
        "title": "a(x,y)",
        "par": 25
    }
    b_data={
        "func": b,
        "title": "b(x,y)",
        "par": 25
    }

    # Para cada una de las funciones graficamos su mapa de calor para comprobar que funcionan
    for funcion_data in [f_data,a_data,b_data]:
        representar_mallado(generar_mallado(funcion_data["func"],parametro=funcion_data["par"]),title=funcion_data["title"])
    plt.show()

if apartado_i:
    
    # En el apartado i buscamos encontrar la solucion a partir de las matrices Tx y Ty, utilizando los parametros predeterminados. Es decir, ejecutar representar_solucion
    representar_solucion(plot_show=True)

if apartado_j:
    # En este apartado comprobaremos lo que sucede al variar C. Para ello, haremos un for loop:
    for numero in [10,20,30]:
        representar_solucion(C=numero,titulo='U(x,y), caso C='+str(numero))
    
    # Hacemos que se generen todos los graficos a la vez
    plt.show()


if apartado_k:
    # Similar al apartado j pero variando kx=ky

    for numero in [0,15,20]:
        representar_solucion(kx=numero,ky=numero, titulo='U(x,y), caso kx=ky='+str(numero))
    
    # Hacemos que se generen todos los graficos a la vez
    plt.show()

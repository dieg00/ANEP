import numpy as np
import matplotlib.pyplot as plt
from figura1 import M as masa_total_analitica


def ajuste_N(N,Neumann=True):
    #Pequeña funcion auxiliar que devuelve N+1 en el caso Neumann y N-1 en Dirichlet
    return N-1+2*Neumann
    

def calcular_matriz(parametro_lambda, parametro_theta,N,Neumann=True,devolver_inversa=True):
    #Esta funcion calcula la matriz (o su inversa), con la que obtendremos los valores de la solucion en el siguiente paso temporal. Requiere los parametros lambda y theta, así como el número de subdivisiones espaciales N. Los toggle booleanos Neumann y devolver_inversa permiten especificar si estamos empleando condiciones de frontera de tipo Neumann, y si queremos obtener la inversa de la matriz; respectivamente.

    #Distinguimos entre los dos casos Neumann y Dirichlet de cara a la forma y valores de la matriz
    multiplicador = 1+Neumann 
    N = ajuste_N(N,Neumann)

    #Computamos los valores de la matriz con 3 diagonales
    matriz = np.diag([1+2*parametro_lambda*parametro_theta]*N)+np.diag([-parametro_theta*parametro_lambda]*(N-1),-1)+np.diag([-parametro_theta*parametro_lambda]*(N-1),1)

    #En el caso Neumann, hay dos entradas que son -2*lambda*theta
    matriz[0,1] *= multiplicador
    matriz[N-1,N-2] *= multiplicador

    #Para comprobar el resultado, se ha añadido una entrada booleana que especifica si se busca la matriz en sí, o su inversa (que al final es lo que usaremos para calcular)
    if devolver_inversa:
        matriz = np.linalg.inv(matriz)
    
    #Junto a otros casos, por problemas de overflow para algunas pruebas se han introducido arrays y variables con datatype float64, que permite mayor precision
    matriz = np.array(matriz,dtype=np.float64)
    return matriz


def generar_mallado(N,M):
    #Genera el mallado de la solución, una matriz de N filas y M columnas. Cada columna representará un instante temporal (las filas serán los distintos valores espaciales). Requiere el número de subdivisiones espaciales N y temporales M.
    mallado = np.zeros([N+1,M+1],dtype=np.float64)
    return mallado


def generar_valores(funcion,x_inicial,x_final,N):
    #Permite popular un array con valores de una función de puntos equiespaciados. Se usa para obtener las condiciones iniciales cuando vienen dadas por una funcion. Requiere la funcion con la que queremos popular el array, los puntos inicial (x_inicial) y final (x_final), y el numero de subdivisiones del intervalo.
    xi = np.linspace(x_inicial,x_final,N+1)
    return funcion(xi)


def generar_mallado_valores(funcion,x_inicial,x_final,N,t_inicial,t_final,M):
    #Como generar_valores pero bidimensional, para el caso del término de creación, en el que se necesita conocer valores en un mallado espaciotemporal. Por tanto, la funcion ha de estar definida f(x,t). Requiere la funcion, los limites de los intervalos espacial y temporal, así como el número de sus subdivisiones.

    #Aprovechamos generar_mallado para crear un array con las dimensiones de la solución u, que compartirá los puntos espaciotemporales en que estamos computando el cálculo numérico.
    f = generar_mallado(N,M)

    #Generamos los arrays espaciales y temporales que formarán la malla
    xi = np.linspace(x_inicial,x_final,N+1)
    ti = np.linspace(t_inicial,t_final,M+1)

    #Recorremos la malla asignando el valor de la funcion en cada punto
    for j in range(M+1):
        for i in range(N+1):
            f[i,j] = funcion(x=xi[i],t=ti[j])
    return f


def inicializar_mallado(condiciones_iniciales,N,M):
    #Esta funcion inicializa el mallado incorporando las condiciones iniciales en la primera columna (correspondiente al primer instante temporal). Requiere el array de condiciones_iniciales, así como las dimensiones del mallado de puntos.

    #Generamos el mallado
    u = generar_mallado(N=N,M=M)

    #Incorporamos las condiciones iniciales
    u[:,0] = condiciones_iniciales
    return u


def calcular_bi(u,i,j,f,N,parametro_lambda,parametro_theta, Delta_t, Neumann = True):
    #Permite computar la entrada i-esima del vector bj, por lo que estamos particularizando un punto de la malla determinado. Requiere la malla u actualizada, las coordenadas i, j; los parametros lambda y theta, el intervalo temporal Delta_t. Además, se añade un toggle booleano entre las dos posibilidades de condiciones de frontera: Neumann y Dirichlet

    #En el caso Neumann, la definicion de bj[0] y bj[N] es distinto al tomar puntos exteriores al mallado que comparten valores con los puntos inmediatamente adyacentes al borde. La unica diferencia en las expresiones es un termino_variable.
    if Neumann and i in [0,N]:

        #Pequeño truco para definir i0=1 si i==0, i0=N-1 si i==N
        i0 = min(i + 1, N - 1)

        #Computa dicho termino_variable
        termino_variable = np.float64(2 * u[i0, j])

    #Si no estamos en el caso Neumann en estos valores en la frontera espacial, tendremos el termino_variable habitual
    else:
        termino_variable = np.float64(u[i+1,j]+u[i-1,j])

    #Computamos la formula con este termino_variable
    return parametro_lambda*(1-parametro_theta)*termino_variable+u[i,j]*(1-2*parametro_lambda*(1-parametro_theta))+Delta_t*(parametro_theta*f[i,j+1]+parametro_theta*f[i,j])

    
def calcular_b(u,j,f,N,parametro_lambda,parametro_theta, Delta_t, Neumann = True):
    #Funcion envoltura de calcular_bi que establece un loop para computar cada entrada del vector b. Requiere la malla u actualizada, la columna temporal j actual; los parametros lambda y theta, el intervalo temporal Delta_t, e incluye un toggle booleano para distinguir en la operativa entre condiciones Neumann y Dirichlet

    #Inicializa el vector b, que medirá N-1 si tenemos condiciones de Dirichlet, y N+1 si tenemos condiciones de Neumann. Nótese que por ello se usa la funcion auxiliar ajuste_N definida anteriormente.
    b = np.zeros([ajuste_N(N,Neumann),1],dtype=np.float64)

    #Itera sobre las distintas entradas a popular del vector i, llamando a calcular_bi en cada caso con el punto espacial correspondiente
    for i in range(N+1):

        #Aqui solo tenemos el problema de definicion de que el vector b tiene entradas de 1 a N-1 en Dirichlet, y de 0 a N en Neumann. Por tanto, queremos NO estar en la situacion Dirichlet AND i=0,i=N
        if Neumann or i not in [0,N]:

            #Si estamos en Dirichlet, tendremos que indexar bj[1] a bj[0] (el array va de 0 a N-2). Para ello usamos el truco de -(not Neumann), que será -1 en este caso, y 0 en Neumann.
            b[i-(not Neumann)]=calcular_bi(u=u,i=i,j=j,f=f,N=N,parametro_lambda=parametro_lambda,parametro_theta=parametro_theta,Neumann=Neumann, Delta_t=Delta_t)

    return b


def iterar_operaciones(u,matriz,N,M,f,parametro_lambda,parametro_theta,Delta_t,Neumann=True):
    #Representa el paso más largo: habiendo configurado las condiciones iniciales y los diversos mallados, es necesario iterar. El proceso será computar b para un tiempo j, y multiplicarle la inversa de la matriz para hallar el estado del sistema en un tiempo j+1.
    #Requeriremos la malla inicializada u, la matriz obtenida con los parametros, los numeros de subdivisiones N y M, el mallado f con los valores del termino de creacion, los parametros lambda y theta, el intervalo Delta_t, y el toggle booleano Neumann para distinguir entre casos.

    #Iteramos para cada tiempo j, que nos permite calcular el sistema en el tiempo j+1 posterior. Por tanto, solo necesitamos iterar hasta j=M-1
    for j in range(M):

        #Se calcula el vector bj
        bj = calcular_b(u=u,j=j,f=f,N=N,parametro_lambda=parametro_lambda,parametro_theta=parametro_theta,Delta_t=Delta_t,Neumann=Neumann)

        #Se multiplica la matriz y el vector para obtener la columna u[:,j+1]. En el caso Dirichlet, tendremos que tener en cuenta que actualizamos solo los valores de la malla entre i==1 e i==N-1. Para ello introducimos "not Neumann:N+Neumann", que es 0:N+1 en Neumann y 1:N en Dirichlet
        u[not Neumann:N+Neumann,j+1]=matriz.dot(bj).ravel()

    return u


def calcular_sistema(parametro_theta,N,M,x_inicial=0,x_final=1,t_inicial=0,t_final=1,Neumann=True,condiciones_son_funcion=True,funcion_condiciones=None, condiciones_iniciales=None, termino_fuente=None):
    #La función que ordena globalmente la lógica de la operación llamando a las distintas funciones anteriores, para generar el mallado de soluciones final u.
    #Requiere el parametro theta, el numero de subdivisiones espaciales N y temporales M, así como los límites de los intervalos espacial ([x_inicial,x_final]) y temporal ([t_inicial,t_final]). El toggle Neumann distingue entre los dos tipos de condiciones de frontera. El toggle condiciones_son_funcion determina si se tomará la funcion funcion_condiciones como funcion para generar el array inicial de condiciones iniciales (en caso de True) o por el contrario se usará el array condiciones_iniciales definidas numericamente (en caso de False). Por ultimo se incluye la funcion de creacion f(x,t) en el termino_fuente
    
    #Se calcula el tamaño de los intervalos espacial y temporal, además del parametro lambda asociado.
    Delta_x = (x_final-x_inicial)/N
    Delta_t = (t_final-t_inicial)/M
    parametro_lambda = Delta_t/(Delta_x**2)

    #Se verifica que se cumple la condicion de convergencia para el lambda y theta resultantes. En caso contrario se genera un error
    if (1-parametro_theta)*parametro_lambda > 1/2:
        print("Yes")
        raise ValueError("Los parametros lambda y theta no cumplen con la condicion de convergencia")
    
    #En caso de que las condiciones iniciales se obtengan a partir de una funcion, se calcula ahora el array condiciones_iniciales a partir de los valores que adopta dicha funcion en los puntos espaciales
    if condiciones_son_funcion:
        condiciones_iniciales=generar_valores(funcion=funcion_condiciones,x_inicial=x_inicial,x_final=x_final,N=N)

    #Se inicializa el mallado con dichas condiciones_iniciales
    u = inicializar_mallado(condiciones_iniciales=condiciones_iniciales,N=N,M=M)

    #Se computa la inversa de la matriz que se utilizará al iterar el método.
    matriz = calcular_matriz(parametro_lambda=parametro_lambda, parametro_theta=parametro_theta,N=N,Neumann=Neumann,devolver_inversa=True)

    #Se computa el mallado con los valores del termino de creacion para cada punto espaciotemporal
    f = generar_mallado_valores(funcion=termino_fuente,x_inicial=x_inicial,x_final=x_final,N=N,M=M,t_inicial=t_inicial,t_final=t_final)

    #Se calcula el estado del sistema en cada tiempo devolviendo la matriz u con toda la evolucion y el estado final
    u = iterar_operaciones(u, matriz, N, M, f, parametro_lambda, parametro_theta, Delta_t, Neumann=Neumann)
    return u


def min_M(parametro_theta,N,t_final,t_inicial=0,x_inicial=0,x_final=1):
    #Esta función computa el minimo valor de M para que pueda converger el sistema. Sigue el criterio de convergencia y requiere N, theta y los limites de los intervalos espacial y temporal.
    return int(np.ceil((2*(1-parametro_theta)*N**2*(t_final-t_inicial))/((x_final-x_inicial)**2)))

##A partir de este punto se encuentra el codigo de nuestro problema en particular, que permite resolver los apartados especificos

def u0(x):
    #Definimos la funcion u0(x) de condiciones iniciales
    return x**2*(1-x)**2

def f(x,t):
    #definimos el termino de creacion
    if t <= np.pi:
        return (1+np.cos(t))*x*(1-x)
    else:
        return 0

def u0_modificada(x):
    #Condiciones iniciales negativas para la ultima pregunta del apartado k
    return -2+x

def f_modificada(x,t):
    #Termino fuente nulo para la ultima pregunta del apartado k
    return 0

def mapacalor(sistema,condiciones):
    #Funcion que genera el mapa de calor, para un sistema, y denotando unas condiciones
    fig, ax = plt.subplots()
    im = ax.imshow(sistema, aspect='auto', cmap='hot')
    im.get_cmap
    cbar = ax.figure.colorbar(im, ax=ax)
    ax.set_title("Sistema u(x,t), condiciones "+condiciones)
    ax.set_ylabel("Eje espacial (x)")
    ax.set_xlabel("Eje temporal (t)")
    y_ticks = np.linspace(0, 1, 11)
    ax.set_yticks(np.linspace(0, N, 11))
    ax.set_yticklabels(np.round(y_ticks,decimals=1))
    x_ticks = np.linspace(0, T, 11)
    ax.set_xticks(np.linspace(0, M, 11))
    ax.set_xticklabels(np.round(x_ticks,decimals=1))
    plt.show()


if __name__ == '__main__':    
    parametro_theta=0.5     #Ajustamos theta, N y T
    N=100
    T=5

    plotear_cosas = True    #Ajusta si queremos que se muestren las figuras, para temas de debugging
    sistemas_normales = True #Ajusta si queremos el calculo de los sistemas normales, para temas de debugging
    sistemas_modificados = True #Ajusta si queremos el calculo de los sistemas modificados, para temas de debugging

    #Ajustamos M. Si usamos un theta distinto de 1, nos guiaremos por el requisito de convergencia. En caso contrario, asignaremos de forma arbitraria M=1000
    if parametro_theta !=1:
        M = min_M(parametro_theta=parametro_theta,N=N,t_final=T)
    else:
        M = 1000

    if sistemas_normales:   #Calculamos el sistema con las condiciones de frontera de Neumann
        print("Running Neumann...")
        sistema_neumann = calcular_sistema(parametro_theta=parametro_theta,N=N,M=M,t_final=T,funcion_condiciones=u0,Neumann=True,termino_fuente=f)
        print("Finished!")

        if plotear_cosas: #Dibujamos un mapa de calor que nos permita ver la evolucion
            mapacalor(sistema_neumann, "Neumann")


            #Calculamos la masa total M(t) integrando u(x,t) respecto a x. Esto es, sumando cada columna y multiplicando por Delta_t
            masa_total_numerica = np.sum(sistema_neumann,axis=0)*(1/N)
            tiempos = np.linspace(0,T,M+1)

            #Representamos la masa total graficamente
            X=np.linspace(0,5,1000)
            Y=[masa_total_analitica(x) for x in X]
            plt.plot(tiempos, masa_total_numerica,c='orange', label='Numérica')
            plt.plot(X,Y, label='Analítica')
            plt.title("Comparativa evolución masa total")
            plt.xlabel("Tiempo t")
            plt.ylabel("Masa total M(t)")
            plt.legend()
            plt.show()


        #Calculamos el sistema con las condiciones de frontera de Dirichlet
        print("Running Dirichlet...")
        sistema_dirichlet = calcular_sistema(parametro_theta=parametro_theta,N=N,M=M,t_final=T,funcion_condiciones=u0,Neumann=False,termino_fuente=f)
        print("Finished!")

        if plotear_cosas:#Dibujamos un mapa de calor que nos permita ver la evolucion
            mapacalor(sistema_dirichlet, "Dirichlet")

        #Obtenemos el valor minimo de la resta para ver si u_N es mayor que u_D
        print("El valor minimo de la resta u_N-u_D es ", np.min(sistema_neumann-sistema_dirichlet))


    if sistemas_modificados:
        #Calculamos los sistemas modificados
        print("Running Neumann modificado...")
        sistema_neumann_modificado = calcular_sistema(parametro_theta=parametro_theta,N=N,M=M,t_final=T,funcion_condiciones=u0_modificada,Neumann=True,termino_fuente=f_modificada)
        print("Now running Dirichlet modificado...")
        sistema_dirichlet_modificado = calcular_sistema(parametro_theta=parametro_theta,N=N,M=M,t_final=T,funcion_condiciones=u0_modificada,Neumann=False,termino_fuente=f_modificada)
        print("Finished!")

        if plotear_cosas:
            mapacalor(sistema_neumann_modificado, "Neumann, modificado")
            mapacalor(sistema_dirichlet_modificado, "Dirichlet, modificado")
        print("El valor minimo de la resta de u_N-u_D (modificados) es ", np.min(sistema_neumann_modificado-sistema_dirichlet_modificado), " y el máximo es ", np.max(sistema_neumann_modificado-sistema_dirichlet_modificado))
        #print("El valor minimo de la resta de u_D-u_N (modificados) es ", np.min(sistema_neumann_modificado-sistema_dirichlet_modificado))



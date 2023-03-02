import numpy as np
import matplotlib.pyplot as plt

ALPHA = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
DELTA = [408, 89, -66, 10, 338, 807, 1238, 1511, 1583, 1462, 1183, 804]
N=len(ALPHA)
c_k = np.fft.rfft(DELTA, norm='forward') #Calculamos la rFFT, que nos devuelve los coeficientes c_k complejos.
print("Los coeficientes c_k son:")
print(c_k)
a_k = 2*c_k.real
b_k = -2*c_k.imag
print(a_k)
print(b_k)
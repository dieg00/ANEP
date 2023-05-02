import numpy as np
import matplotlib.pyplot as plt
A = np.zeros([3,6])
B = np.sum(A,axis=0)
print(A.shape)
print(np.linspace(1,3.5,5))

m=1
delta_t=0.2323
plt.plot([1,2],[3,4])
plt.text(1.5, 3.5, f"Tiempo: {m * delta_t:.2f}", fontsize=22)
plt.show()
import numpy as np
import matplotlib.pyplot as plt
X=np.linspace(0,1.5*np.pi,1000)
def M(x):
    if x<=np.pi:
        return 1/30+1/6*(x+np.sin(x))
    else:
        return 1/30+np.pi/6
Y=[M(x) for x in X]
Y_pi = np.linspace(np.min(Y)+0.05,M(np.pi),1000)
X_pi = np.pi*np.ones_like(Y_pi)
plt.plot(X,Y)
plt.scatter(np.pi,M(np.pi),c='black')
plt.xlabel("Tiempo t")
plt.ylabel("Masa total M(t)")
plt.title("Evolución masa total")
plt.plot(X_pi, Y_pi, linestyle='dashed', c='black')
plt.text(np.pi, np.min(Y), "t=π", ha='center', va='bottom', fontsize=12)
plt.show()
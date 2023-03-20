import numpy as np
import matplotlib.pyplot as plt
X=np.linspace(0,1.5*np.pi,1000)
def M(x):
    if x<=np.pi:
        return 1/30+1/6*(x+np.sin(x))
    else:
        return 1/30+np.pi/6
Y=[M(x) for x in X]
plt.plot(X,Y)
plt.scatter(np.pi,M(np.pi))
plt.show()
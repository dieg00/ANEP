import numpy as np
print(np.zeros(2))


a = [1,2,3,4,5]



output = np.diag(a[:-1])+np.diag(a[1:])-np.diag(a[1:-1],-1)-np.diag(a[1:-1],1)
print(output)
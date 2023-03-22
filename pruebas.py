import numpy as np

A = np.zeros([3,6])
B = np.sum(A,axis=0)
print(B.shape)
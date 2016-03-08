import numpy as np
from numpy.linalg import inv


matrix = np.array([
    [0.28591, 0.00001, 0.00000, 0.00000],
    [0.00000, 0.91794, 0.00577, 0.00000],
    [0.00303, 0.00913, 0.79342, 0.00405],
    [0.00102, 0.00000, 0.01198, 0.99471],


])

err_matrix = np.array([
    [0.0029,  0.0000,  0.0000,  0.0000],
    [0.0000,  0.0009,  0.0003,  0.0000],
    [0.0004,  0.0003,  0.0014,  0.0002],
    [0.0002,  0.0000,  0.0004,  0.0002],

])

print("Matrix")
print(matrix)
inv_matrix = inv(matrix)
print("inverse matrix")
print(inv_matrix)

print("Error Matrix")
print(err_matrix)
inv_err_mat = - np.dot(inv_matrix, np.dot(err_matrix, inv_matrix))
print(inv_err_mat)


"""
crazy shit coming up here
"""
e = np.zeros((7, 4))
e[0] = [20,93,97,2500]
e[1] = [81,282,266,6650]
e[2] = [91,389,345,8770]
e[3] = [631,2879,2398,66863]
e[4] = [114,577,484,12995]
e[5] = [37,263,251,6257]
e[6] = [56,292,255,6982]
print("True events")
print(inv_matrix.dot(e.T).T)
print("error")
print(inv_err_mat.dot(e.T).T)



import numpy as np
from numpy.linalg import inv
from util import np_array_to_latex

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
print(np_array_to_latex(matrix,2))
inv_matrix = inv(matrix)
print("inverse matrix")
print(inv_matrix)
print(np_array_to_latex(inv_matrix,2))


print("Error Matrix")
print(err_matrix)
inv_err_mat = - np.dot(inv_matrix, np.dot(err_matrix, inv_matrix))
print(inv_err_mat)
print(np_array_to_latex(inv_err_mat,2))


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
print('data from experiment')
print(np_array_to_latex(e,0))

print("True events")
ntrue = inv_matrix.dot(e.T).T
print(inv_matrix.dot(e.T).T)
print(np_array_to_latex(inv_matrix.dot(e.T).T,3))
print("error")
dntrue = inv_err_mat.dot(e.T).T
print(dntrue)
print(np_array_to_latex(inv_err_mat.dot(e.T).T,3))

fixes = np.array([
    [0.09]*3+[2.0],
    [0.2]*3+[4.3],
    [0.36]*3+[7.7],
    [0.52]*3+[10.8],
    [0.22]*3+[4.7],
    [-0.01]*3+[-0.2],
    [-0.08]*3+[-1.6],
], dtype=np.float32)

print('fixes_latex')
print(np_array_to_latex(fixes,2))
print('L')
lumi_errs = [[4.249604],[5.691792],[4.454466],[16.43293],[4.848926],[4.276552],[6.104764]]
lumi = [[463.979], [667.5236], [486.7642], [2246.568], [535.908], [450.6], [709.698]]
print(np_array_to_latex(lumi,4))
print('ERROR')
print(np_array_to_latex(lumi_errs))
lumis = np.array([a*4 for a in lumi], dtype=np.float32)
sigma = e/lumis+fixes
print('sigma')
print(sigma)
print(np_array_to_latex(sigma,5))

dsigma = np.empty_like(sigma)
for i in range(0,7):
    for j in range(0,4):
        dsigma[i][j] = np.sqrt((1.0/lumi[i][0]*dntrue[i][j])**2+(-ntrue[i][j]/(lumi[i][0])**2*lumi_errs[i][0])**2)

print('dsigma')
print(np_array_to_latex(dsigma,4))

"""
differential cross-setion
"""







import numpy as np
from scipy.linalg import solve
import random as rand
import matplotlib.pyplot as plt
import numpy.linalg as npl

def density_h (Z, r):
    #Modified for task 2
    ns = np.exp(-Z*r)/np.pi
    return ns

Z = 2
N  = 200
A  = np.zeros((N,N))
ns = np.zeros((N,1))

length = 10
r, h= np.linspace(0,length,N+2, retstep=True)
"""
A[0][0] = 0
A[0,0] = 0
A[0][1] = -r[0] / r[1] /(1-Z*h)
A[N-1][N-1] = 0
for i in range(1, N-1):
    A[i][i-1] =-0.5/h**2
    A[i][i]   =1.0/h**2-1.0/r[i]
    A[i][i+1] =-0.5/h**2
"""
A[0,0] = 1.0/h**2-1.0/r[1]
A[0,1] = -0.5/h**2

A[N-1,N-1] = 1.0/h**2-1.0/r[N]
A[N-1,N-2] = -0.5/h**2
for i in range(1, N-1):
    A[i][i-1] =-0.5/h**2
    A[i][i]   =1.0/h**2-1.0/r[i+1]
    A[i][i+1] =-0.5/h**2


x,vec = npl.eig(A)
y = np.argmin(x)
eigen_vec = vec[:,y]
if (eigen_vec[10] < 0):
    eigen_vec = -1.0*eigen_vec
scale = np.trapz(eigen_vec*eigen_vec, r[2::])
ev_norm = eigen_vec / np.sqrt(scale)
print('Analytical value', -13.6/27.211) # 27.211 EV = 1HT, 1 HT 
print('Smallest eigenvalue = ', x[y])
#print(vec)
print(ev_norm.shape)
print(ev_norm[1])
print(ev_norm[0])
ev_norm = np.insert(ev_norm, 0,666)
print(ev_norm.shape)
print(ev_norm[1])
print(ev_norm[0])
ev_norm = np.append(ev_norm, 667)
print(ev_norm.shape)
print(ev_norm[-1])
vh_real = []
"""
for i in range(1,N):
    #vh_real.append(1.0/r[i] -(1.0+1.0/r[i])*np.exp(-2.0*r[i]))
    vh_real.append(2.0* np.exp(-r[i]))
    #vh_simulated.append(x[i]/r[i])
plt.plot(r[3::], vh_real, label='Real $V_h$')
plt.plot(r,ev_norm/r, 'gx', label='Corresponding eigenvector / r')
plt.xlabel('r [a.u]')
plt.ylabel('Potential [a.u]')
plt.legend()
plt.show()
"""

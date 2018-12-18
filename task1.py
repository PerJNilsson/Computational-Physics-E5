import numpy as np
from scipy.linalg import solve
import random as rand
import matplotlib.pyplot as plt

def density_h (Z, r):
    ns = np.exp(-2*Z*r)/np.pi
    return ns


N  = 100
A  = np.zeros((N,N))
ns = np.zeros((N,1))

length = 20
r, h= np.linspace(0,length,N, retstep=True)

A[0][0] = 1
A[N-1][N-1] = 1
for i in range(1, N-1):
    A[i][i-1] =1/h**2
    A[i][i]   =-2/h**2
    A[i][i+1] =1/h**2
A[1][0] = 0
# Get the electron density with BC
Z =1
for i in range(1,N-1):
    ns[i] = -4*np.pi*r[i]*density_h(Z,r[i])
ns[0] = 0
ns[N-1] = 1

x = solve(A,ns)

vh_real = []
vh_simulated = []

for i in range(1,N):
    vh_real.append(1.0/r[i] -(1.0+1.0/r[i])*np.exp(-2.0*r[i]))
    
    vh_simulated.append(x[i]/r[i])

plt.plot(r[1::], vh_real, label='Real $V_h$')
plt.plot(r[1::], vh_simulated, 'g*', label='Solved x')
plt.xlabel('r [a.u]')
plt.ylabel('Potential [a.u]')
plt.legend()
plt.show()


import numpy as np
import random as rand
import matplotlib.pyplot as plt
import numpy.linalg as npl

def density_h (Z, r):
    ns = np.exp(-2*Z*r)/np.pi
    return ns

old_E = 10000.0
energy_groundstate = 0
N  = 1600
A  = np.zeros((N,N))
ns = np.zeros((N,1))
ground_energy = 0
length = 5
r, h= np.linspace(0.01,length,N, retstep=True)
print(h)
# Initial guess
# Get the electron density with BC
Z =2
for i in range(1,N-1):
    ns[i] = 100#-4*np.pi*r[i]*density_h(Z,r[i])
ns[0] = 0
ns[N-1] = 1

while (np.abs(old_E - energy_groundstate) > 1E-8):
    A[0][0] = 1
    A[N-1][N-1] = 1
    for i in range(1, N-1):
        A[i][i-1] =1/h**2
        A[i][i]   =-2/h**2
        A[i][i+1] =1/h**2


    x = npl.solve(A,ns)

    A[0,0] = 0
    A[0][1] = -r[0] / r[1] /(1-Z*h)
    A[N-1][N-1] = 1
    for i in range(1, N-1):
        A[i][i-1] = -0.5/h**2
        A[i][i] = 1.0/h**2 - 2.0/r[i] + x[i]/r[i]
        A[i][i+1] = -0.5/h**2


    
    eig, vec = npl.eig(A)
    y = np.argmin(eig)
    eigen_vec = vec[:,y]
    if (eigen_vec[10] < 0): # in case it negative
        eigen_vec = 1.0*eigen_vec
    scale = np.trapz(eigen_vec*eigen_vec, r) #integrating with trapzoid method
    ground_energy = eig[y]
    ev_norm = eigen_vec / np.sqrt(scale)
    wave_fun = 1.0 / (np.sqrt(4*np.pi))*ev_norm/r
    charge_density = wave_fun*wave_fun

    # Update density
    for i in range(1,N-1):
        ns[i] = -4*np.pi*r[i]*charge_density[i]
    ns[0] = 0
    ns[-1] = -1/np.sqrt(4*np.pi)

    old_E = energy_groundstate
    tmp_vec = np.zeros(N)
    for i in range(0,N):
        tmp_vec[i] = x[i]*charge_density[i] *r[i] *4*np.pi

    energy_groundstate = 2*ground_energy -np.trapz(tmp_vec,r)


tmp_vec = np.zeros(N)
for i in range(0,N):
    tmp_vec[i] = x[i]*charge_density[i]/r[i]

energy_groundstate = 2*ground_energy -np.trapz(tmp_vec,r)
print('Energy Hartree: ',energy_groundstate, 'a.u.')
print('CFA Z=2:         -2.750 a.u.')
print('CFA Z=27/16:     -2.848 a.u.')
plt.plot(r, 4*np.pi*r*r*charge_density, 'g--', mfc='none')
plt.plot(r, 4*r*r*(2)**3*np.exp(-2*2*r), label='Unscreened Z=2')
plt.plot(r, 4*r*r*(27/16)**3*np.exp(-2*(27/16)*r), label='Screened Z=27/16')
plt.xlabel('r / [a.u.]')
plt.ylabel(r'$\rho(r)$')
plt.title('Probability density functions')
plt.show()

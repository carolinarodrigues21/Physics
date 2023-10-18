# Partícula Livre em 1D - Oscilador Harmônico
# Carolina Niklaus - 14/6/2022

import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import * 
from scipy import special

k = 1
m = 10000

def V(k,x):
    return k * x**2 /2

def niveisDeEnergia(omega, n):
    En = []
    for i in n:
        E = hbar * omega * (i + 1/2) * 10**(36)
        En.append(E)
    return En

def x0(m,omega):
    return np.sqrt(hbar * 10**(36)/(m*omega))

def estadosEstacionarios(x,n,x_o):
    Hn = special.hermite(n)
    psi = 1/(np.pi**(1/4)*x_o**(1/2)) * 1/(np.sqrt(2**n * special.factorial(n))) * Hn(x/x_o) * np.exp(-x**2/(2*x_o**2))
    return psi


x = np.linspace(-4,4,100)
omega = np.sqrt(k/m)

n = np.arange(0,7)
x_o = x0(m,omega)
niveis = niveisDeEnergia(omega,n)

estado0 = estadosEstacionarios(x,0,x_o)
estado1 = estadosEstacionarios(x,1,x_o)
estado2 = estadosEstacionarios(x,2,x_o)
estado3 = estadosEstacionarios(x,3,x_o)


# print("------Diferença entre os níveis de energia------------")
# d = 0
# while (d<6):
#     print(niveis[d+1] - niveis[d])
#     d+=1
# print("-----------------------------------------------------")

#Gráfico Potencial com níveis de energia

# plt.plot(x,V(k,x))
# plt.xlabel("x")
# plt.ylabel("V(x)")
# # l = 0
# # for j in niveis:
# #     cor = ['red','blue', 'green', 'pink', 'purple', 'orange', 'black']
# #     plt.axhline(y=j, c=cor[l], label = "E(%i) = %5f" %(l,niveis[l]))
# #     l+=1
# # plt.plot(x, estados[0])
# plt.plot(x,estado0)
# plt.plot(x,estado1)
# plt.plot(x,estado2)
# plt.plot(x,estado3)
# plt.grid(True)
# # plt.legend()
# plt.show()


fig, axs = plt.subplots(2, 1)
axs[0].plot(x,V(k,x))
axs[0].set_xlabel("x")
axs[0].set_ylabel("V(x)")
l = 0
for j in niveis:
    cor = ['red','blue', 'green', 'pink', 'purple', 'orange', 'black']
    axs[0].axhline(y=j, c=cor[l], label = "E(%i) = %5f" %(l,niveis[l]))
    l+=1
axs[0].grid(True)
axs[0].legend(loc='upper right')

axs[1].plot(x,estado0, label = r'$\psi_0$', color = 'red')
axs[1].plot(x,estado1, label = r'$\psi_1$', color = 'blue')
axs[1].plot(x,estado2, label = r'$\psi_2$', color = 'green')
axs[1].plot(x,estado3, label = r'$\psi_3$', color = 'pink')
axs[1].grid(True)
axs[1].legend(loc='upper right')
plt.show()
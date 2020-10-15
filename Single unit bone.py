import numpy as np
import matplotlib.pyplot as plt

# variable decleration 
dens_end = ((0.26+0.75*3)/4)*1000
c = 500
dt = 0.001
t = 100
N = np.int(t/dt)

f = int(input("Enter force in N "))
a = int(input("Input area in m^2 "))
p_start = int(input("Enter start density in kg/m^3 "))
p = np.zeros(N)
p[0]=p_start


#function decleration
def strainenergy(sigma,p):
    e_strain=sigma**2/(6850*p**1.49)  #sigma*epsilon
    return e_strain

def deltap(u,k):
    dp = c*(u-k)*dt
    return dp
    
def sigma_function(f,a):
    sigma = f/a
    return sigma


def pcalc(f,p,dt,t,N):
    sigma = sigma_function(f,a)
    k = strainenergy(sigma,dens_end)
    print(k)
    for i in range(N-1):
        u = strainenergy(sigma,p[i])
        p_dot = deltap(u,k)
        p[i+1]=p[i]+p_dot
    print(u)
    return p

pcalc(f,p,dt,t,N)
tlin = np.linspace(0,t,N)
plt.close('all')
plt.figure()
plt.plot(tlin,p,label='$density$',linestyle='-')
plt.xlabel('$days$')
plt.ylabel('$density [kg/m^3]$')
plt.title('Density graph with start density ' + str(p_start) + 'kg/m^3, force of ' + str(f) + 'N and area ' + str(a)+ 'm^2') 
plt.grid()
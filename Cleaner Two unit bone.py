import numpy as np
import matplotlib.pyplot as plt

# variable decleration 
dens_end = ((0.26*1+0.75*4)/4)*1000#apparent density between 0.26 and 0.75 gg/cm^3 according to literature for femoral neck bone
c = 25000 #this was a random guess, affects the rate of change drastically however
dt = 0.001
t = 100
N = np.int(t/dt)
sigma = [0,0] #reserving memory for both sigma, k and pdt
k = [0,0]
pdt = [0,0]
u_k_info = [[],[]]

f_in1=[]
f_in2=[]

#function decleration
def forcedis(f,a1,a2,p1,p2):
    f1 = (f*(p1/(p1+p2))+f*(a1/(a1+a2)))/2 #figured the force was both dependent on the density and the area of the force
    f2 = (f*(p2/(p1+p2))+f*(a2/(a1+a2)))/2
    return [f1,f2]

def strainenergy(sigma,p):
    e_strain=sigma**2/(6850*p**1.49) #e = sigma*epsilon, epsilon = sigma/E,  E = sigma^2/epsilon=stress^2/strain
    #E (denuminator) defined as femoral neck bone according to literature Trabecular bone modulus–density relationships depend on anatomic site
    return e_strain

def deltap(u,k):
    dp = c*(u-k)*dt #change in density depends on current strain energy (u), goal strain energy (k) and rate on bone change (c)
    return dp
    
def sigma_function(f,a):
    sigma = f/a #by definition of strain
    return sigma

def pcalc(a1, a2,f,p_1,p_2,dt,t,N):
    #fk = forcedis(f,a[0],a[1],p_1[0],p_2[0])
    for i in range(N-1): #for each time step, we are going to re-evaluate the density
        f_dis=forcedis(f,a[0],a[1],p_1[i],p_2[i]) #the force on each element is dependent on density
        for j in range(2):
            sigma[j]=sigma_function(f_dis[j],a[j]) #the current stress on each element
            u_1=strainenergy(sigma[0],p_1[i]) #strain energy element 1 at iteration i
            u_2=strainenergy(sigma[1],p_2[i]) #strain energy element 2 at iteration i
            pdt[0]=deltap(u_1,k) #change in density for element 1
            pdt[1]=deltap(u_2,k) #change in density for element 2
        p_1[i+1]=p_1[i]+pdt[0] #new density of element 1 for index i+1
        p_2[i+1]=p_2[i]+pdt[1] #new density of element 2 for index i+1
    return p_1,p_2
 
#the actual code
while(1):
    
    #user inputs of force, area of element 1 and 2 and start densities
    f = int(input("Enter force in N "))
    a = [float(input("Input area of bone element 1 in m^2 "))]
    a.append(float(input("Input area of bone element 2 in m^2 ")))
    p_start = [float(input("Enter start density of element 1 in kg/m^3 "))]
    p_start.append(float(input("Enter start density of element 2 in kg/m^3 ")))
    k = strainenergy(0.5*f, dens_end)#reference strainenergy K
    
    #preparing data and reserving memory
    p_1 = np.zeros(N)
    p_1[0]=p_start[0]
    
    p_2 = np.zeros(N)
    p_2[0]=p_start[1]
    
    #calculation
    pcalc(a[0],a[1],f,p_1,p_2,dt,t,N)
    #standard matplot stuff, nothing too scary or original
    tlin = np.linspace(0,t,N)
    plt.close('all')
    plt.figure()
    plt.plot(tlin,p_1,label='density element 1',linestyle='-')
    plt.plot(tlin,p_2,label='density element 2', linestyle = '-')
    plt.legend()
    plt.xlabel('$days$')
    plt.ylabel('$density [kg/m^3]$')
    plt.title('Density graph with a force of ' +str(f)+ 'N, \n start density element 1 ' + str(p_start[0]) + 'kg/m^3 and area ' + str(a[0])+ 'm^2 \nstart density element 2 ' + str(p_start[1]) + 'kg/m^3 and area ' + str(a[1])+ 'm^2 ') 
    plt.grid()
    plt.show()
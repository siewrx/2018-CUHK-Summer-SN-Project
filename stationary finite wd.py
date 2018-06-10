import numpy as np
from mpmath import mp,mpf
mp.dps=80
import matplotlib.pyplot as plt
# define the constants
h,G,pi,m_p,m_e,c=(mpf('6.62607004e-34'),mpf('6.67408e-11'),mp.pi,mpf('1.672621898e-27'),mpf('9.10938356e-31'),mpf('299792458'))
k=6*G*m_p/(m_e*c**2)
k1=pi*(m_e**4)*(c**5)/(3*h**3)
k2=(h/(c*m_e))*mp.cbrt(3/(16*pi*m_p))
M_s=mpf('1.989e30')
r_s=695700000
rho_c=5*10**13
rho=[rho_c]
P=[]
dr=200
r=[0.0]
M=[0.0]
R=[]
M1=[]
def rad(r):
    r.append(r[-1]+dr)

def mass(M,rho,r):
    if len(M)==1: 
        M.append(mpf('4/3')*pi*rho[-1]*(dr**3))
    else:
        # since density decreases outward, rho is overestimated. Hence the radius should be underestimated.
        M.append(M[-1]+4*dr*pi*rho[-1]*r[-2]**2)       

def dens(rho,M,r):
    F=k2*mp.cbrt(rho[-1])
    rho.append(rho[-1]-dr*k*((M[-1]+M[-2])/2)*rho[-1]*mp.sqrt(1+F**2)/(r[-1]*F)**2)

def pres(P,rho):
        F=k2*mp.cbrt(rho[-1])
        P.append(k1*np.real(F*mp.sqrt(1+F**2)*(-3+2*F**2)+3*mp.log(F+mp.sqrt(1+F**2))))  
# produce the central pressure
pres(P,rho)
#start producing data
def wd(rho,P,r,M):
    while rho[-1]>=0:
        rad(r)
        mass(M,rho,r)
        dens(rho,M,r)
        pres(P,rho)
    # Mass-radius curve
    plt.plot(r,M)
    plt.xlabel('R (m)')
    plt.ylabel('M (kg)') 
    plt.title('Mass profile')
    plt.show()
    # Pressure-radius curve
    plt.plot(r,P)
    plt.xlabel('R (m)')
    plt.ylabel('P (Pa)') 
    plt.title('Pressure profile')
    plt.show()
    # Rho-radius curve
    plt.plot(r,rho)
    plt.xlabel('R (m)')
    plt.ylabel(r'$\rho\ (kg/m^3)$')
    plt.title('Density profile') 
    plt.show()

    
def wd1(rho,P,r,M):
    while rho[-1]>=0 and P[-1]>=0:
        rad(r)
        mass(M,rho,r)
        dens(rho,M,r)
        pres(P,rho)
def rM(n,rhostep,R,M1,rho,r,M,P):
    while len(R)<=n:
        rho=[rho[0]+rhostep]
        P=[k*rho[0]**(4/3)]
        r=[0.0]
        M=[0.0]  
        wd1(rho,P,r,M)
        R.append(r[-1]/r_s)
        M1.append(M[-1]/M_s)
    plt.plot(M1,R)
    plt.xlabel(r'$Mass\ (M_s)$')
    plt.ylabel(r'$radius\ (R_s)$') 
    plt.title('Radius vs Mass')
    plt.show()
    
import time
start = time.time()
wd(rho,P,r,M)
#R.append(r[-1]/r_s)
#M1.append(M[-1]/M_s)  
#rM(50,5*10**9,R,M1,rho,r,M,P)
end = time.time()
print(end - start)  
print('Radius='+str(r[-1]/r_s)+' '+'solar radius')
print('Mass='+str(M[-1]/M_s)+' '+'solar mass')   



   
   


#Problem 4.13
'''   Develop a program that applies the fourth-order Runge–Kutta and bisection methods to solve
    the eigenvalue problem of thestationary one-dimensional Schrodinger equation.
        Find the two lowest eigenvalues and eigenfunctions for an electron in the potential well '''

import math
import matplotlib.pyplot as plt
import numpy as np
import cmath

#lovest energy
x_0=5
hbar=1
m=1
L=x_0
V_0=5
k=3*(math.pi)/L

E_0=2*k**2/2

def F(x,y,p):
    k2=((2*m)/(hbar**2))*(E_0-V(x,x_0,V_0))
    return k2*y

def runge_kutta_4_second_order(F,x,y,p,x_final,h=0.05):
    p_data=[p]
    y_data=[y]
    x_data=[x]
    while x<x_final:
        j1=h*F(x,y,p)
        k1=h*p
        j2=h*F(x+h/2,y+k1/2,p+j1/2)
        k2=h*(p+j1/2)
        j3=h*F(x+h/2,y+k2/2,p+j2/2)
        k3=h*(p+j2/2)
        j4=h*F(x+h,y+k3,p+j3)
        k4=h*(p+j3)
        p=p+(1/6)*(j1+2*j2+2*j3+j4)
        p_data.append(p)
        y=y+(1/6)*(k1+2*k2+2*k3+k4)
        y_data.append(y)
        x=x+h
        x_data.append(x)
    return x_data,y_data,p_data

def k(x):
    return cmath.sqrt( (2*m)/(hbar**2)*(E_0- V(x,x_0,V_0))).real


def V(x,x_0,V_0):
    ''' This function is given potential condition in the problem set.'''
    if x > 0 and x < x_0:
        return V_0*(x/x_0)
    else:
        return V_0

x_potential=np.linspace(-5,10,100) # generating x data for potential values
V_potential=[]
for i in x_potential:
    ''' Generating potential y data'''
    V_potential.append(V(i,x_0,V_0))


p=(math.pi/L)/15
x,y,p=runge_kutta_4_second_order(F,-5,0,p,10)
p2=(math.pi/L)/30
x2,y2,p2=runge_kutta_4_second_order(F,-5,0,p2,10)

fig, ax = plt.subplots()
plt.xlabel("x domain",fontsize=12)
plt.ylabel("U(x) and V(x)",fontsize=12)

#Graphicing by matplotlib.pyplot
ax.grid(True, which='both')
ax.axhline(y=0, color='k')
ax.axvline(x=0, color='k')

plt.plot(x_potential,V_potential,color='blue')# plotting potential.
plt.plot(x,y,color="red")
plt.plot(x2,y2,color='green')
plt.show()
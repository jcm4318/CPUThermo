# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 17:13:08 2020

@author: james
"""

#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


width = 20 # width of microprocessor in mm
height = 2 # height of microprocessor in mm
h = 0.5
hx = h # x step size
hy = h # y step size
Ta = 0
initval = 100 # initial guess for T at all points
q = 0.5
k = 150


#%% Plotting the domain

x = np.arange(-width/2 - hx, width/2 + 2*hx, hx)
y = np.arange(-height/2 - hx, height/2 + 2*hy, hy)
X, Y = np.meshgrid(x, y)
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)

fig1 = plt.figure(figsize = (width,height))
plt.plot(X,Y,'.',color = 'red',ms = 5)

# for i in range(len(x)):
#     for j in range(len(y)):
#         plt.text(x[i]+(x[1]-x[0])/8,y[j]+(y[1]-y[0])/8,"(%i, %d)" %( i , j ))

plt.title('Microprocessor Discrete Points')
plt.tight_layout()
plt.show()


#%%

T = np.zeros((len(y),len(x)))
T.fill(initval)

updatemat = np.zeros((len(y),len(x)))

Q = np.zeros((len(y),len(x)))
Q.fill(q)

def nextij(T,Q,i,j,h = h,k = k):
    # print(X[i,j],Y[i,j])
    # print(T[i,j])
    # print()
    
    # now update the relevant ghost point(s) using the boundary conditions:
    if (i == 1): # left side
        Hc = 1.31 * (abs(T[i,j]) - Ta)**(1/3)
        updatemat[i-1,j] = (h*Hc/k) * (T[i,j] - Ta) + T[i,j+1] - T[i,j-1] + T[i+1,j]

    if (j == 1): # top
        Hc = 1.31 * (abs(T[i,j]) - Ta)**(1/3)
        updatemat[i,j-1] = (h*Hc/k) * (T[i,j] - Ta) + T[i,j+1] - T[i-1,j] + T[i+1,j]
    
    if (i == T.shape[0] - 2): # right
        Hc = 1.31 * (abs(T[i,j]) - Ta)**(1/3)
        updatemat[i+1,j] = (h*Hc/k) * (Ta - T[i,j]) - T[i-1,j] + T[i,j-1] - T[i,j+1]
        
    if (j == T.shape[1] - 2): # bottom
        Hc = 1.31 * (abs(T[i,j]) - Ta)**(1/3)
        updatemat[i,j+1] = (h*Hc/k) * (Ta - T[i,j]) - T[i,j+1] + T[i,j-1] - T[i+1,j]
    
    # then return the value for the point in question
    return (1/4) * ( (h**2/k) * Q[i,j] - ( T[i-1,j] + T[i,j-1] + T[i+1,j] + T[i,j+1] ) )
    
def nextiter():
    for i in range(1,T.shape[0]-1): # for every point in the domain of T
        for j in range(1,T.shape[1]-1):
            # print("i,j: ")
            # print(i)
            # print(j)
            # print()
            updatemat[i,j] = nextij(T,Q,i,j) # calculate and assign the T value to updatemat
            # print(updatemat[i,j])
            
    for i in range(0,T.shape[0]): # for every point in the domain
        for j in range(0,T.shape[1]):
            T[i,j] = updatemat[i,j] # move it from updatemat to T
            

for i in range(50): # run the code for 200 iterations
    nextiter()
    
#%%
#Tt = np.transpose(T)

fig = plt.figure(figsize = (width,height))
plt.title('Microprocessor Discrete Points')
# plt.tight_layout()
# plt.plot(X,Y,'.',color = 'red',ms = 5)
# for i in range(len(x)):
#     for j in range(len(y)):
#         plt.text(x[i]+(x[1]-x[0])/8,y[j]+(y[1]-y[0])/8, T[i,j] )

#fig = plt.figure()
ax = fig.gca(projection='3d')
# Plot the surface.

Xplot = X[ 1:X.shape[0]-2 , 1:X.shape[1]-2 ] # trim out the ghost points
Yplot = Y[ 1:Y.shape[0]-2 , 1:Y.shape[1]-2 ]
Tplot = T[ 1:T.shape[0]-2 , 1:T.shape[1]-2 ]


surf = ax.plot_surface(Xplot, Yplot, Tplot, cmap=cm.binary,
                        linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-110, 110)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
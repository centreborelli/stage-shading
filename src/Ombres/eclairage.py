import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pylab as plt
from scipy.linalg import norm
import copy as c

#os.chdir('C:\\Users\\dugas\\Desktop\\stage')
#A = np.load('fuji.npy')
plt.figure()
plt.clf()

X,Y = np.meshgrid(np.linspace(0,1,31),np.linspace(0,1,31))
tableaubase = ThreeBumps(X,Y)
# os.chdir('C:\\Users\\dugas\\Documents\\stage-shading\\data\\heights')
# tableaubase = np.load("deltaplane.npy")
n,m = np.shape(tableaubase)
# A = np.zeros((3*n,3*m))
# A[n:2*n,m:2*m] = tableaubase
A = tableaubase

def SimuleOmbre(height,params) :

    n,m = np.shape(height)
    alpha,beta,gamma,h = params

    if alpha != 0 :
        sgn1 = -int(alpha/np.abs(alpha))
    if beta != 0 :
        sgn2 = -int(beta/np.abs(beta))

    if beta != 0 :
        sb = int(beta/abs(beta))
    else :sb = 1

    if alpha != 0 :
        sa = int(alpha/abs(alpha))
    else : sa = 1

    # print(sgn1,sgn2)
    O = []
    Ombrage = np.ones((n,m))
    Traite = np.zeros((n,m))

    for i in range(-min(0,-sb*(n-1)),-max(-n,-sb*(-1)),-sb)  :
        for j in range(-min(0,-sa*(n-1)),-max(-n,-sa*(-1)),-sa) :

            if not Traite[i,j] :

                if alpha != 0 :

                    f = lambda t : beta/alpha *(t-i)+j
                    for k in range(i+sgn1,max(sgn1*n,sgn1),sgn1) :
                        if int(f(k)) >= 0 :
                            try  :
                                if (height[k,int(f(k))] < height[i,j] - np.abs((k-i)*h/alpha*gamma)) and np.abs(int(f(k)) - f(k)) < 1/2 :
                                    Ombrage[k,int(f(k))] = 0
                                    Traite[k,int(f(k))] = 1
                            except IndexError :
                                pass
                            try :
                                if (height[k,int(f(k))+1] < height[i,j] - np.abs((k-i)*h/alpha*gamma)) and np.abs(int(f(k))+1 - f(k)) < 1/2 :
                                    Ombrage[k,int(f(k))+1] = 0
                                    Traite[k,int(f(k))+1] = 1
                            except IndexError :
                                pass


                if beta != 0 :

                    f = lambda y : alpha/beta  * (y-j) + i
                    for k in range(j+sgn2,max(sgn2*n,sgn2),sgn2) :
                        if int(f(k)) >= 0 :
                            try :
                                if (height[int(f(k)),k] < height[i,j] - np.abs((k-j)*h/beta*gamma)):# and np.abs(int(f(k)) - f(k)) < 1/2 :
                                    Ombrage[int(f(k)),k] = 0
                                    Traite[int(f(k)),k] = 1
                            except IndexError :
                                pass
                            try :
                                if (height[int(f(k))+1,k] < height[i,j] - np.abs((k-j)*h/beta*gamma)):# and np.abs(int(f(k))+1 - f(k)) < 1/2 :
                                    Ombrage[int(f(k))+1,k] = 0
                                    Traite[int(f(k))+1,k] = 1
                            except IndexError :
                                pass
            x = c.deepcopy(Ombrage)
            O.append(x)

    return(Ombrage,O)



def GenerateRHS(height,params,ombre = False):
    α,β,γ,h = params
    hx,hy = np.gradient(height,h)
    Intensity = (-α*hx-β*hy+γ)/np.sqrt(1+hx**2+hy**2)
    Omega = height>0
    Ombre = 0
    if ombre:
        Ombre,O = SimuleOmbre(height,params)
        Intensity = Intensity #* SimuleOmbre(height,params)
    return Intensity*Ombre,Omega,Ombre,O

alpha = 5
beta= 0
gamma = 1
alpha, beta, gamma = alpha/norm([alpha,beta,gamma],2) , beta/norm([alpha,beta,gamma],2),gamma/norm([alpha,beta,gamma],2)

tableauinter,oui,ombre,O = GenerateRHS(A,(alpha,beta,gamma,1/30),ombre = True)


fig,ax = plt.subplots(subplot_kw = {"projection":"3d"})
surf = ax.plot_surface(X,Y,A,cmap = cm.coolwarm,linewidth = 0,antialiased = False)
fig.colorbar(surf,shrink = 0.5,aspect = 5)

plt.show()

#
# plt.imshow(tableauinter,cmap = 'gray')
# plt.show()

#
# plt.figure()
# plt.clf()
# for i in range(len(O)) :
#     if i %5000 == 0 :
#         print(i)
#         plt.imshow(O[i])
#         plt.pause(0.1)
# plt.show()


# C = 1
# #T = [t/10*np.pi for t in range(11)]
# T = [1]
# Itot = []
# gamma = 1
#
# for theta in T :
#
#     alpha = np.array([np.cos(gamma)*np.sin(theta),np.cos(theta),-np.sin(gamma)*np.sin(theta)])
#
#     x = np.shape(A)[0]
#     y = np.shape(A)[1]
#
#     dx= 1
#     dy = 1
#
#     Px = np.zeros((x,y,1))
#     Py = np.zeros((x,y,1))
#
#     Px[:-1,:] = (A[1:,:] - A[:-1,:])/dx
#     Px[-1,:] = (A[-1,:] - A[-2,:])/dx
#
#     Py[:,:-1] = (A[:,1:] - A[:,:-1])/dy
#     Py[:,-1] = (A[:,-1] - A[:,-2])/dy
#
#     u = np.sqrt(1+Px[:,:,0]**2+Py[:,:,0]**2)
#
#     N = [[(-Px[i,j,0]/u[i,j],-Py[i,j,0]/u[i,j],1/u[i,j]) for j in range(y)] for i in range(x)]
#
#     N = np.array(N)
#
#     I = np.zeros((x,y))
#
#     for i in range(x) :
#         for j in range(y) :
#             I[i,j] = C*np.sum(N[i,j]*alpha)
#     Itot.append(I)
#
#
# for I in Itot :
#     plt.imshow(I)
#     plt.pause(0.1)
# plt.show()
#
# tableauinter = Itot[0]
import numpy as np
import matplotlib.pyplot as plt

def Terrier(x,y) :
    r2 = (x-0.5)**2+(y-0.5)**2
    r2max=0.09
    danscercle = (r2<r2max)
    return (y>0.5)*np.sqrt((1-r2/r2max)*danscercle)
#
# h=1/100
# X,Y = np.meshgrid(np.arange(0,1,h),np.arange(0,1,h),indexing='ij')
# height = Terrier(X,Y)
# masque = np.zeros((100,100))
# masque[21:80,50] = 2
# masque[21:80,51] = 4
# plt.figure("Relief originale")
# plt.clf()
# plt.imshow(height,cmap=plt.cm.RdBu_r)
# plt.colorbar()
# plt.figure("Masque")
# plt.clf()
# plt.imshow(masque,cmap=plt.cm.RdBu_r)
# plt.colorbar()
def GenerateRHS_discontinue(height,masque,params):
    α,β,γ,h = params
    hx,hy = np.gradient(height,h)
    (n,m) = height.shape
    for i in range(n) :
        for j in range(m) :
            if masque[i,j] == 1 :
                hy[i,j] = (height[i+1,j]-height[i,j])/h
            elif masque[i,j] == 2 :
                hx[i,j] = (height[i,j]-height[i,j-1])/h
            elif masque[i,j] == 3 :
                hy[i,j] = (height[i,j]-height[i-1,j])/h
            elif masque[i,j] == 4 :
                hx[i,j] = (height[i,j+1]-height[i,j])/h
    Intensity = (α*hx+β*hy+γ)/np.sqrt(1+hx**2+hy**2)
    return Intensity


def Pont(x,y) :
    a=8
    p = 1-a*(x-0.5)**2
    return (p>0)*(y>=0.4)*(y<=0.6)*p

h=1/100
X,Y = np.meshgrid(np.arange(0,1,h),np.arange(0,1,h),indexing='ij')
height = Pont(X,Y)
masque = np.zeros((100,100))
masque[15:86,39] = 2
masque[15:86,40] = 4
masque[15:86,60] = 2
masque[15:86,61] = 4



Etats = np.ones((n,m)) # 0 = Far ; 1 = Trial ; 2 = Accepted
Etats[1:-1,1:-1] = np.zeros((n-2,m-2))
Bord = Etats == 1

U0 = np.zeros((n,m))
U0[:,:] = height
U0[1:-1,1:-1] = np.inf*np.ones((n-2,m-2))
V = np.zeros((n,m)) # Valeurs des trials
V[:,:] = U0
α,β,γ=0,0,1
params = α,β,γ,h

I = GenerateRHS_discontinue(height,masque,params)
# I2,Omega = GenerateRHS(height,params)
# [U] = FMM_discontinue(Etats,V,U0,I,Bord,params,fc.solve_quad,masque)
# [U2] = FMM_discontinue(Etats,V,U0,I2,Bord,params,fc.solve_quad,masque)

plt.show()
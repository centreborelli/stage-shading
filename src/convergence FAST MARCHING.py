import numpy as np
import pylab as plt
import os
path = 'D:/L3/Stage/modules python'
os.chdir(path)
import fonctions_FMM as fm
import fonctions_calculs as fc
import fonctions_espace as fe

plt.ion()
plt.show()

H = []
Err = []
α,β,γ= fe.parametres(0,90)


for k in range(10,16) :
    print(k)
    h=(1.4)**(-k)
    H.append(h)
    params = α,β,γ,h
    X,Y = np.meshgrid(np.arange(0,1,h),np.arange(0,1,h),indexing='ij')
    height = fe.OneBump(X,Y)
    (n,m) = height.shape
    CAS = np.zeros((n,m))
    # height=np.zeros((n,m))
    """On impose zéro aux bords du rectangle, et on résoud SFS avec le fast marching en partant du bord"""

    Etats = np.ones((n,m)) # 0 = Far ; 1 = Trial ; 2 = Accepted
    Etats[1:-1,1:-1] = np.zeros((n-2,m-2))
    Bord = Etats == 1

    U = np.zeros((n,m)) # valeur de la fonction
    U[:,:] = height
    U[1:-1,1:-1] = np.inf*np.ones((n-2,m-2))
    V = np.zeros((n,m)) # Valeurs des trials
    V[1:-1,1:-1] = np.inf*np.ones((n-2,m-2))

    # On peut aussi imposer zéro aux bords du relief et commencer la FMM sur ce bord en utililsant la ligne suivante
    # U,V,Etats = Initialise_contour(height,0)

    X = np.arange(0, n, 1)
    Y = np.arange(0, m, 1)
    X, Y = np.meshgrid(X, Y)
    I,Omega = fe.GenerateRHS(height,params)

    ## FAST MARCHING

    [U,CAS,ITERATION] = fm.FFM(Etats, V, U, I, Bord, params, fc.solve_quad, Cas=True,iterations=True)

    Err.append(np.sum(np.abs(U-height))*h**2)


plt.figure(1)
plt.clf()
[A,B] = np.polyfit(np.log10(H),np.log10(Err),1)
Reg = H**A*10**B
plt.loglog(H,Err,'o')
plt.loglog(H,Reg)
print(A)
plt.show()






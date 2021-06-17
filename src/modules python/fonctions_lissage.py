import numpy as np

# import pylab as plt
import os
path = 'C:/Users/dugas/Documents/stage-shading/src/modules python'
os.chdir(path)
import fonctions_FMM as fm
import fonctions_calculs as fc
import fonctions_espace as fe
import fonctions_lissage as fl
import copy as c
import scipy.spatial as scp
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.linalg import norm
from scipy.linalg import solve

def lisse_moy(V) :
    A = c.deepcopy(V)
    n,m = V.shape
    for i in range(n) :
        for j in range(m) :
            try :
                A[i,j] = np.mean(np.array([V[i-1,j],V[i,j],V[i+1,j],V[i-1,j-1],V[i,j-1],V[i+1,j-1],V[i-1,j+1],V[i,j+1],V[i+1,j+1]]))
                # A[i,j] = (1/9)*(V[i-1,j]+V[i,j]+V[i+1,j]+V[i-1,j-1]+V[i,j-1]+V[i+1,j-1]+V[i-1,j+1]+V[i,j+1]+V[i+1,j+1])
            except IndexError :
                pass
    return(A)

def lisse_med(V) :
    A = c.deepcopy(V)
    n,m = V.shape
    for i in range(n) :
        for j in range(m) :
            try :
                A[i,j] = np.median(np.array([V[i-1,j],V[i,j],V[i+1,j],V[i-1,j-1],V[i,j-1],V[i+1,j-1],V[i-1,j+1],V[i,j+1],V[i+1,j+1]]))
            except IndexError :
                pass
    return(A)

def env_conv(n,m,X,Y,V) :

    eps = 10**(-3)
    D = np.zeros((n*m,3))
    for i in range(n*m) :
        D[i,:] = np.array([X[i//m,i%m],Y[i//m,i%m],V[i//m,i%m]])
    A = scp.ConvexHull(D)

    B = np.zeros((n,m))

    for i in range(n) :
        for j in range(m) :
            u = np.array([j,i,1])
            for simplex in A.simplices :
                Mat = np.array([D[simplex,0],D[simplex,1],[1,1,1]])
                try :
                    a = solve(Mat,u)
                    z = float(n-1)

                    k = Mat - np.array([[z,  0.,  0.],[z, z,  0.],[ 1.,  1.,  1.]])
                    l = Mat - np.array([[z, z,  0.],[z,  0.,  0.],[ 1.,  1.,  1.]])
                    p = Mat - np.array([[z,  0.,  z],[z, z,  0.],[ 1.,  1.,  1.]])
                    q = Mat - np.array([[0, z,  0.],[z,  0.,  0.],[ 1.,  1.,  1.]])

                    if (a[0] >= -eps) and (a[1] >=-eps) and (a[2] >= -eps) and (k.nonzero()[0].size >0 or k.nonzero()[1].size > 0) and (l.nonzero()[0].size >0 or l.nonzero()[1].size > 0) and (p.nonzero()[0].size >0 or p.nonzero()[1].size > 0) and (q.nonzero()[0].size > 0 or q.nonzero()[1].size > 0) :

                        B[i,j] = a[0]*D[simplex[0],2] + a[1]*D[simplex[1],2] + a[2]*D[simplex[2],2]
                        break

                except np.linalg.LinAlgError :
                    # print('ok')
                    None
            if i == 30 and j == 40 :
                None

    #
    # plt.figure()
    # plt.clf()
    #
    # fig,ax = plt.subplots(subplot_kw = {"projection":"3d"})
    # surf = ax.plot_surface(X,Y,B,cmap = cm.coolwarm,linewidth = 0,antialiased = False)
    # fig.colorbar(surf,shrink = 0.5,aspect = 5)
    #
    # plt.show()

    return(B)


def lisse_conv(V) :
    n,m = V.shape
    W = np.zeros((n,m))

    for i in range(0,n,5) :
        for j in range(0,m,5) :
            try:
                K = V[i-2:i+3,j-2:j+3]
                X0 = np.arange(0, 5, 1)
                Y0 = np.arange(0, 5, 1)
                X0, Y0 = np.meshgrid(X0, Y0)
                W[i-2:i+3,j-2:j+3] = fl.env_conv(5,5,X0,Y0,K)
            except (IndexError,scp.qhull.QhullError) :
                W[i,j] = V[i,j]
    return(W)


def gaussian2d(m,n,std):
    return [[1 / (std**2 * 2*np.pi) * np.exp(-float(x - n//2)**2/(2*std**2)) * np.exp(-float(y - n//2)**2/(2*std**2)) for y in range(n)] for x in range(m)]

def lisse_conv_gauss(V) :

    M,N = V.shape
    std = 2
    g = gaussian2d(M,N, std)
    gg = np.fft.fftshift(g)
    smooth = np.real(np.fft.ifft2(np.fft.fft2(gg) * np.fft.fft2(V)))

    return(smooth)

def dev2(V,i,j) :
    return(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1]-4*V[i,j])


def liss_dev2(V) :
    n,m = V.shape
    W = np.zeros((n,m))
    for i in range(n):
        for j in range(m) :
            try :
                u = (1/9) * (np.sum(np.array([[dev2(V,k,l) for k in range(i-1,i+2)] for l in range(j-1,j+2)])))
                W[i,j] = 1/4*(V[i+1,j]+V[i-1,j]+V[i,j+1]+V[i,j-1] - u)
            except IndexError :
                None

    return(W)
























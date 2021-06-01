import numpy as np
import os
path = 'D:/L3/Stage/modules python'
os.chdir(path)
import fonctions_calculs as fc
import pylab as plt


def voisins(i,j,n,m) :
  "renvoie [[voisins_de_gauche, voisins_de_droite],[voisins_du_haut, voisins_de_bas]]"
  if i == 0 :
      if j == 0 :
          return([[(0,1)],[(1,0)]])
      elif j==m-1 :
          return([[(0,m-2)],[(1,m-1)]])
      else :
        return([[(0,j-1),(0,j+1)],[(1,j)]])
  elif i == n-1 :
      if j == 0 :
          return([[(n-1,1)],[(n-2,0)]])
      elif j==m-1 :
          return([[(n-1,m-2)],[(n-2,m-1)]])
      else :
          return([[(n-1,j-1),(n-1,j+1)],[(n-2,j)]])
  else :
      if j == 0 :
          return([[(i,1)],[(i-1,0),(i+1,0)]])
      elif j==m-1 :
          return([[(i,m-2)],[(i-1,m-1), (i+1,m-1)]])
      else :
          return([[(i,j-1),(i,j+1)],[(i-1,j), (i+1,j)]])

def voisins_per(i,j,n,m) :
  "renvoie [[voisins_de_gauche, voisins_de_droite],[voisins_du_haut, voisins_de_bas]]"
  if i == 0 :
      if j == 0 :
          return([[(0,m-1),(0,1)],[(n-1,0),(1,0)]])
      elif j==m-1 :
          return([[(0,m-2),(0,0)],[(n-1,m-1),(1,m-1)]])
      else :
        return([[(0,j-1),(0,j+1)],[(n-1,j),(1,j)]])
  elif i == n-1 :
      if j == 0 :
          return([[(n-1,m-1),(n-1,1)],[(n-2,0),(0,0)]])
      elif j==m-1 :
          return([[(n-1,m-2),(n-1,0)],[(n-2,m-1),(0,m-1)]])
      else :
          return([[(n-1,j-1),(n-1,j+1)],[(n-2,j),(0,j)]])
  else :
      if j == 0 :
          return([[(i,m-1),(i,1)],[(i-1,0),(i+1,0)]])
      elif j==m-1 :
          return([[(i,m-2),(i,0)],[(i-1,m-1), (i+1,m-1)]])
      else :
          return([[(i,j-1),(i,j+1)],[(i-1,j), (i+1,j)]])

def voisins_gradient(U,pi,pj) :
    '''Renvoie les coordonnées des voisins à choisir pour le calcul du nouveau point, en indiquant s'ils ont été pris à gauche ou à droite (en haut ou en bas) avec epx et epsy'''
    (n,m) = U.shape
    res = voisins(pi,pj,n,m)
    if len(res[0]) == 1 :
        q_gradx = res[0][0]
        if q_gradx[0] < pi :
            epsx = -1
        else :
            epsx = 1
    else :
        vois_gauche, vois_droite = res[0][0],res[0][1]
        if U[vois_gauche[0],vois_gauche[1]] < U[vois_droite[0],vois_droite[1]] :
            q_gradx = vois_gauche
            epsx = -1 # calculer la dérivée dans le bon sens dans solve_quad
        else :
            q_gradx = vois_droite
            epsx = 1
    if len(res[1]) == 1 :
        q_grady = res[1][0]
        if q_grady[1] < pj :
            epsy = -1
        else :
            epsy = 1
    else :
        vois_haut, vois_bas = res[1][0],res[1][1]
        if U[vois_haut[0],vois_haut[1]] < U[vois_bas[0],vois_bas[1]] :
            q_grady = vois_haut
            epsy = -1
        else :
            q_grady = vois_bas
            epsy = 1
    return(q_gradx, q_grady, epsx, epsy)


def voisins_gradient_per(U,pi,pj) :
    '''Renvoie les coordonnées des voisins à choisir pour le calcul du nouveau point, en indiquant s'ils ont été pris à gauche ou à droite (en haut ou en bas) avec epx et epsy'''
    (n,m) = U.shape
    res = voisins_per(pi,pj,n,m)

    vois_gauche, vois_droite = res[0][0],res[0][1]
    if U[vois_gauche[0],vois_gauche[1]] < U[vois_droite[0],vois_droite[1]] :
        q_gradx = vois_gauche
        epsx = -1 # calculer la dérivée dans le bon sens dans solve_quad
    else :
        q_gradx = vois_droite
        epsx = 1

    vois_haut, vois_bas = res[1][0],res[1][1]
    if U[vois_haut[0],vois_haut[1]] < U[vois_bas[0],vois_bas[1]] :
        q_grady = vois_haut
        epsy = -1
    else :
        q_grady = vois_bas
        epsy = 1
    return(q_gradx, q_grady, epsx, epsy)



def Initialise_contour(height,seuil) :
    '''renvoie les tableaux initiaux pour une FMM qui commence sur les contours d'un relief
    Fonctionne pour des formes strictement incluses dans le rectangle'''
    V = np.inf*np.ones_like(height)
    U = np.inf*np.ones_like(height)
    M = 1*(height > seuil)
    Mx,My = np.gradient(M)
    B = np.hypot(Mx,My)>0
    Etats = 1*B
    U[B] = 0
    V[B] = 0
    return(U,V,Etats)

def FFM_var(Etats,V,U,I,Bord,params,solve_quad,variation,cas = False, iterations = False,plot = []) :
    (n,m) = U.shape
    Z = 2*np.ones((n,m))
    it = 0
    ITERATION = np.zeros((n,m))
    CAS = np.zeros((n,m))
    if len(plot)>0 :
        plt.ion()
        plt.show()
        from matplotlib import cm
    while not((Etats==Z).all()) :
        k = np.argmin(V)
        i = k//m
        j = k - i*m # (i,j) sont les coordonnées du point Trial de valeur minimale
        Etats[i,j] = 2
        if 1 in plot :
            plt.figure("plot Etats")
            plt.clf()
            plt.imshow(Etats,cmap=plt.cm.RdBu_r)
            plt.colorbar()
            plt.pause(0.01)
        if 2 in plot :
            plt.figure("plot U")
            plt.clf()
            plt.imshow(U,cmap=plt.cm.RdBu_r)
            plt.colorbar()
            plt.pause(0.01)
        ITERATION[i,j] = it
        it+=1
        V[i,j] = np.inf
        vois = voisins(i,j,n,m)
        list_voisins = vois[0]+vois[1]
        for p in list_voisins :
            pi,pj = p
            if Etats[pi,pj] != 2 and not(Bord[pi,pj]) :
                Etats[pi,pj] = 1
                q_gradx, q_grady, epsx, epsy = voisins_gradient_per(U,pi,pj)
                (r1,r2),cas = solve_quad(U[q_gradx[0],q_gradx[1]],U[q_grady[0],q_grady[1]],I[pi,pj],epsx,epsy,params)
                CAS[pi,pj] = cas
                if variation(U,r1,i,j,pi,pj) < variation(U,r2,i,j,pi,pj) :
                    u = r1
                else :
                    u = r2
                U[pi,pj] = u # actualisation de la valeur de p
                V[pi,pj] = u
    res = [U]
    if cas :
        res.append(CAS)
    if iterations :
        res.append(ITERATION)
    return(res)


def FFM(Etats,V,U,I,Bord,params,solve_quad,cas = False, iterations = False,plot = []) :
    (n,m) = U.shape
    Z = 2*np.ones((n,m))
    it = 0
    ITERATION = np.zeros((n,m))
    CAS = np.zeros((n,m))
    if len(plot)>0 :
        plt.ion()
        plt.show()
        from matplotlib import cm
    while not((Etats==Z).all()) :
        k = np.argmin(V)
        i = k//m
        j = k - i*m # (i,j) sont les coordonnées du point Trial de valeur minimale
        Etats[i,j] = 2
        if 1 in plot :
            plt.figure("plot Etats")
            plt.clf()
            plt.imshow(Etats,cmap=plt.cm.RdBu_r)
            plt.colorbar()
            plt.pause(0.01)
        if 2 in plot :
            plt.figure("plot U")
            plt.clf()
            plt.imshow(U,cmap=plt.cm.RdBu_r)
            plt.colorbar()
            plt.pause(0.01)
        ITERATION[i,j] = it
        it+=1
        V[i,j] = np.inf
        vois = voisins(i,j,n,m)
        list_voisins = vois[0]+vois[1]
        for p in list_voisins :
            pi,pj = p
            if Etats[pi,pj] != 2 and not(Bord[pi,pj]) :
                Etats[pi,pj] = 1
                q_gradx, q_grady, epsx, epsy = voisins_gradient_per(U,pi,pj)
                u,cas = solve_quad(U[q_gradx[0],q_gradx[1]],U[q_grady[0],q_grady[1]],I[pi,pj],epsx,epsy,params)
                CAS[pi,pj] = cas
                U[pi,pj] = u # actualisation de la valeur de p
                V[pi,pj] = u
    res = [U]
    if cas :
        res.append(CAS)
    if iterations :
        res.append(ITERATION)
    return(res)
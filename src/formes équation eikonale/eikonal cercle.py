import numpy as np
import os
path = 'D:/L3/Stage/Fast marching SFS/CERCLE'
os.chdir(path)




def voisins(i,j) :
  "renvoie [[voisins_de_gauche, voisins_de_droite],[voisins_du_haut],[voisins_de_bas]]"
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



h=1


import numpy as np
import pylab as plt

plt.ion()
plt.show()

n = 201
m = 201
i0,j0 = 100, 100 # centre du carré
r = 80

def est_dans_cercle(i,j,r) :
    return((i-i0)**2+(j-j0)**2 < r**2)
def est_au_bord(i,j,r) :
    L = voisins(i,j)
    for voisin in L[0] :
        if est_dans_cercle(voisin[0],voisin[1],r) :
            return(True)
    for voisin in L[1] :
        if est_dans_cercle(voisin[0],voisin[1],r) :
            return(True)
    return(False)

Etats = 2*np.ones((n,m)) ## n = m ## 0 = Far ; 1 = Trial ; 2 = Accepted
U = np.zeros((n,m))
V = np.inf*np.ones((n,m)) # Valeurs des trials



for i in range(n) :
    for j in range(n) :
        if est_dans_cercle(i,j,r) :
            Etats[i,j] = 0
            U[i,j] = np.inf # valeur de la fonction
        elif est_au_bord(i,j,r) :
            Etats[i,j] = 1
            U[i,j] = 0
            V[i,j] = 0


# TODO : en démarrant sur un pixel random, on n'arrive pas à calculer tous les points, mais en initialisant le contour du bord à 1 ça marche


np.save("cercle Etats initial",Etats)
np.save("cercle V initial",V)
np.save("cercle U initial",U)


Z = 2*np.ones((n,m))




def solve_quad(x,y) :
    if x == np.inf :
        return(y+h) # on considère qu'on a un maximum sur l'axe des x
    elif y == np.inf :
        return(x+h)
    else :
        b = -2*(x+y)
        c = x**2+y**2-h**2
        delta = b**2-8*c
        return((-b+np.sqrt(delta))/4)

def voisins_gradient(pi,pj) :
    res = voisins(pi,pj)
    if len(res[0]) == 1 :
        q_gradx = res[0][0]
    else :
        vois_gauche, vois_droite = res[0][0],res[0][1]
        if U[vois_gauche[0],vois_gauche[1]] < U[vois_droite[0],vois_droite[1]] :
            q_gradx = vois_gauche
        else :
            q_gradx = vois_droite
            # print(U[vois_droite[0],vois_droite[1]],U[vois_gauche[0],vois_gauche[1]])
    if len(res[1]) == 1 :
        q_grady = res[1][0]
    else :
        vois_haut, vois_bas = res[1][0],res[1][1]
        if U[vois_haut[0],vois_haut[1]] < U[vois_bas[0],vois_bas[1]] :
            q_grady = vois_haut
        else :
            q_grady = vois_bas
    return(q_gradx, q_grady)



## Algo :

plt.figure(2)

while not((Etats==Z).all()) :
    k = np.argmin(V)
    i = k//m
    j = k - i*m # (i,j) sont les coordonnées du point Trial de valeur minimale
    Etats[i,j] = 2
    # plt.clf()
    # plt.imshow(Etats)
    # plt.pause(0.01)
    V[i,j] = np.inf
    vois = voisins(i,j)
    list_voisins = vois[0]+vois[1]
    for p in list_voisins :
        pi,pj = p
        if Etats[pi,pj] != 2 :
            Etats[pi,pj] = 1
            q_gradx, q_grady = voisins_gradient(pi,pj)
            u = solve_quad(U[q_gradx[0],q_gradx[1]], U[q_grady[0],q_grady[1]])
            U[pi,pj] = u # actualisation de la valeur de p
            V[pi,pj] = u


# np.save("cercle eikonal",U)
plt.clf()
plt.imshow(U)
plt.colorbar()


plt.show()
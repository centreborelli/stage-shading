import numpy as np
import os
path = 'D:/L3/Stage/Fast marching SFS/DELTAPLANE'
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

n = 101
m = 101
i0,j0 = 50, 50 # centre du carré
a = 20
b = 35
i1,j1 = i0,j0+a
i2,j2 = i0-b,j0-2*a
i3,j3 = i0,j0-a
i4,j4 = i0+b,j0-2*a


def est_dans_deltaplane(i,j) :
    if i<= i0 :
        return((j-j2)*(i1-i2)<=(j1-j2)*(i-i2) and (j-j2)*(i3-i2)>=(j3-j2)*(i-i2))
    else :
        return((j-j1)*(i4-i1)<=(j4-j1)*(i-i1) and (j-j3)*(i4-i3)>=(j4-j3)*(i-i3))


def est_au_bord(i,j) :
    L = voisins(i,j)
    for voisin in L[0] :
        if est_dans_deltaplane(voisin[0],voisin[1]) :
            return(True)
    for voisin in L[1] :
        if est_dans_deltaplane(voisin[0],voisin[1]) :
            return(True)
    return(False)

Etats = 2*np.ones((n,m)) ## n = m ## 0 = Far ; 1 = Trial ; 2 = Accepted
U = np.zeros((n,m))
V = np.inf*np.ones((n,m)) # Valeurs des trials


for i in range(n) :
    for j in range(m) :
        if est_dans_deltaplane(i,j) :
            Etats[i,j] = 0
            U[i,j] = np.inf # valeur de la fonction
        elif est_au_bord(i,j) :
            Etats[i,j] = 1
            U[i,j] = 0
            V[i,j] = 0

Bord = Etats ==1


# np.save("deltaplane Etats initial",Etats)
# np.save("deltaplane V initial",V)
# np.save("deltaplane U initial",U)

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

plt.figure(1)
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
        if Etats[pi,pj] != 2 and not(Bord[pi,pj]) :
            Etats[pi,pj] = 1
            q_gradx, q_grady = voisins_gradient(pi,pj)
            u = solve_quad(U[q_gradx[0],q_gradx[1]], U[q_grady[0],q_grady[1]])
            U[pi,pj] = u # actualisation de la valeur de p
            V[pi,pj] = u




# np.save("deltaplane eikonal",U)

plt.clf()
plt.imshow(U,cmap=plt.cm.RdBu_r)
plt.colorbar()

plt.show()
import numpy as np
import os
path = 'D:/L3/Stage/Fast marching SFS/RECTANGLE'
os.chdir(path)

h=1


import numpy as np
import pylab as plt

plt.ion()
plt.show()

n = 200
m = 300
Etats = np.ones((n,m)) # 0 = Far ; 1 = Trial ; 2 = Accepted
Etats[1:-1,1:-1] = np.zeros((n-2,m-2))

U = np.zeros((n,m)) # valeur de la fonction
U[1:-1,1:-1] = np.inf*np.ones((n-2,m-2))
V = np.inf*np.ones((n,m)) # Valeurs des trials


np.save("rectangle Etats initial",Etats)
np.save("rectangle V initial",V)
np.save("rectangle U initial",U)

Z = 2*np.ones((n,m))

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


plt.figure(1)
plt.clf()

## Algo :

debut = True
while not((Etats==Z).all()) :
    k = np.argmin(V)
    i = k//m
    j = k - i*m # (i,j) sont les coordonnées du point Trial de valeur minimale
    Etats[i,j] = 2
    # plt.clf()
    # plt.imshow(U)
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


# np.save("rectangle eikonal",U)
plt.imshow(U)
plt.colorbar()

plt.show()
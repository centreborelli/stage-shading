import numpy as np
import pylab as plt


plt.ion()
plt.show()



## Fonctions annexes


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





RESIDUS = [] # stock les parties imaginaires des racines des polynomes du second degré

def solve_quad(x, y,cp,epsx,epsy) :
    '''Renvoie la racine la plus grande de l'équation du second degré à résoudre, dépendant du nombre de voisin disponibles'''
    if x == np.inf :
        a = cp**2-β**2
        b = 2*(-cp**2*y+β**2*y+h*β*γ*epsy)
        c = (cp**2*(h**2+y**2)-h**2*γ**2-2*h*γ*β*y*epsy-β**2*y**2)
    elif y == np.inf :
        a = cp**2-α**2
        b = 2*(-cp**2*x+α**2*x+h*α*γ*epsx)
        c = (cp**2*(h**2+x**2)-h**2*γ**2-2*h*γ*α*x*epsx-α**2*x**2)
    else :
        a = 2*cp**2-(α*epsx+β*epsy)**2
        b = 2*(-cp**2*(x+y) + (α*epsx+β*epsy)*(α*x*epsx+β*y*epsy+h*γ))
        c = cp**2*(h**2+x**2+y**2)-(γ*h+α*x*epsx+β*y*epsy)**2
    delta = b**2-4*a*c
    RESIDUS.append((np.roots([a,b,c]).imag))
    return(np.max((np.roots([a,b,c]).real)))



def voisins_gradient(pi,pj) :
    '''Renvoie les coordonnées des voisins à choisir pour le calcul du nouveau point, en indiquant s'ils ont été pris à gauche ou à droite (en haut ou en bas) avec epx et epsy'''
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

## Algo :

## Initialisation en chargeant relief

# TODO : penser à prendre h identique à celui de l'image enregistré
h=1 # c'est toujours h=1 pour les reliefs obtenus par résolution de l'équation eikonale

height = np.load("../data/heights/cône.npy")
(n,m) = height.shape

"""On impose zéro aux bords du rectangle, et on résoud SFS avec le fast marching en partant du bord"""

Etats = np.ones((n,m)) # 0 = Far ; 1 = Trial ; 2 = Accepted
Etats[1:-1,1:-1] = np.ones((n-2,m-2))
U = np.zeros((n,m)) # valeur de la fonction
U[1:-1,1:-1] = np.inf*np.ones((n-2,m-2))
V = np.zeros((n,m)) # Valeurs des trials
V[1:-1,1:-1] = np.inf*np.ones((n-2,m-2))


X = np.arange(0, n, 1)
Y = np.arange(0, m, 1)
X, Y = np.meshgrid(X, Y)



## Choix de la lumière et récupération de l'intensité (lignes issues du Jupyter)

# Shading parameters
ω = 0.10*np.array([1,2]) # c'est sur le réel devant le tableau qu'il faut jouer pour obtenir des images convainquantes (il détermine le poids des composantes en x et y  par rapport à celle en z).
α,β,γ = np.append(ω,1)/np.linalg.norm(np.append(ω,1))

def GenerateRHS(height,params):
    α,β,γ,h = params
    hx,hy = np.gradient(height,h)
    Intensity = (α*hx+β*hy+γ)/np.sqrt(1+hx**2+hy**2)
    Omega = height>0
    return Intensity,Omega

params = α,β,γ,h
I,Omega = GenerateRHS(height,params)





## FAST MARCHING

(n,m) = height.shape
Z = 2*np.ones((n,m))

while not((Etats==Z).all()) :
    k = np.argmin(V)
    i = k//m
    j = k - i*m # (i,j) sont les coordonnées du point Trial de valeur minimale
    Etats[i,j] = 2
    V[i,j] = np.inf
    vois = voisins(i,j,n,m)
    list_voisins = vois[0]+vois[1]
    for p in list_voisins :
        pi,pj = p
        if Etats[pi,pj] != 2 :
            Etats[pi,pj] = 1
            q_gradx, q_grady, epsx, epsy = voisins_gradient(pi,pj)
            u = solve_quad(U[q_gradx[0],q_gradx[1]],U[q_grady[0],q_grady[1]],I[pi,pj],epsx,epsy)
            U[pi,pj] = u # actualisation de la valeur de p
            V[pi,pj] = u


plt.figure(1)
plt.clf()
plt.imshow(height)
plt.colorbar()
plt.figure(2)
plt.clf()
plt.imshow(I,cmap='gray')
plt.colorbar()
plt.figure(3)
plt.clf()
plt.imshow(U)
plt.colorbar()
I2,Omega = GenerateRHS(U,params)
plt.figure(4)
plt.clf()
plt.imshow(I2,cmap='gray')
plt.colorbar()
plt.show()

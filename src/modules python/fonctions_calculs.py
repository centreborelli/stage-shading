import numpy as np
import os
path = 'D:/L3/Stage/modules python'
os.chdir(path)
import fonctions_FMM as fm

def solve_quad(x, y, cp,epsx,epsy,params) :
    '''Renvoie la racine la plus grande de l'équation du second degré à résoudre, dépendant du nombre de voisin disponibles, en indiquant comment elle a été calculée'''
    α,β,γ,h = params
    if x == np.inf :
        a = cp**2-β**2
        b = 2*(-cp**2*y+β**2*y+h*β*γ*epsy)
        c = (cp**2*(h**2+y**2)-h**2*γ**2-2*h*γ*β*y*epsy-β**2*y**2)
        cas = 1
    elif y == np.inf :
        a = cp**2-α**2
        b = 2*(-cp**2*x+α**2*x+h*α*γ*epsx)
        c = (cp**2*(h**2+x**2)-h**2*γ**2-2*h*γ*α*x*epsx-α**2*x**2)
        cas = 2
    else :
        a = 2*cp**2-(α*epsx+β*epsy)**2
        b = 2*(-cp**2*(x+y) + (α*epsx+β*epsy)*(α*x*epsx+β*y*epsy+h*γ))
        c = cp**2*(h**2+x**2+y**2)-(γ*h+α*x*epsx+β*y*epsy)**2
        cas = 3
    delta = b**2-4*a*c
    return((np.max(np.roots([a,b,c])).real),cas,delta)



def solve_quad2(x, y, cp,epsx,epsy,params) :
    '''Renvoie les racines de l'équation du second degré à résoudre, dépendant du nombre de voisin disponibles, en indiquant comment elle a été calculée'''
    α,β,γ,h = params
    if x == np.inf :
        a = cp**2-β**2
        b = 2*(-cp**2*y+β**2*y+h*β*γ*epsy)
        c = (cp**2*(h**2+y**2)-h**2*γ**2-2*h*γ*β*y*epsy-β**2*y**2)
        cas = 1
    elif y == np.inf :
        a = cp**2-α**2
        b = 2*(-cp**2*x+α**2*x+h*α*γ*epsx)
        c = (cp**2*(h**2+x**2)-h**2*γ**2-2*h*γ*α*x*epsx-α**2*x**2)
        cas = 2
    else :
        a = 2*cp**2-(α*epsx+β*epsy)**2
        b = 2*(-cp**2*(x+y) + (α*epsx+β*epsy)*(α*x*epsx+β*y*epsy+h*γ))
        c = cp**2*(h**2+x**2+y**2)-(γ*h+α*x*epsx+β*y*epsy)**2
        cas = 3
    delta = b**2-4*a*c
    return((np.roots([a,b,c]).real),delta)

def variation_per(U,r,i,j,pi,pj) :
    (n,m) = U.shape
    [[voisins_g, voisins_d],[voisins_h, voisins_b]]=fm.voisins_per(i,j,n,m)
    if i==pi :
        if pj-j == 1 :
            if U[voisins_g[0],voisins_g[1]] == np.inf :
                return(0)
            else :
                return(abs(U[voisins_g[0],voisins_g[1]]-2*U[i,j]+U[voisins_d[0],voisins_d[1]]))
        else : #pj-j=-1
            if U[voisins_d[0],voisins_d[1]] == np.inf :
                return(0)
            else :
                return(abs(U[voisins_g[0],voisins_g[1]]-2*U[i,j]+U[voisins_d[0],voisins_d[1]]))
    else :
        if pi-i == 1 :
            if U[voisins_h[0],voisins_h[1]] == np.inf :
                return(0)
            else :
                return(abs(U[voisins_h[0],voisins_h[1]]-2*U[i,j]+U[voisins_b[0],voisins_b[1]]))
        else : #pi-i=-1
            if U[voisins_b[0],voisins_b[1]] == np.inf :
                return(0)
            else :
                return(abs(U[voisins_h[0],voisins_h[1]]-2*U[i,j]+U[voisins_b[0],voisins_b[1]]))

def variation_non(U,r,i,j,pi,pj) :
    return(0)
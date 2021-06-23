import numpy as np


def OneBump(x,y):
    bump = 0.5-3.*((x-0.5)**2+(y-0.5)**2)
    return np.maximum(bump, np.zeros_like(x))
def ThreeBumps(x,y):
    bump1 = 0.3-3*((x-0.4)**2+(y-0.5)**2)
    bump2 = 0.25-3*((x-0.65)**2+(y-0.6)**2)
    bump3 = 0.25-3*((x-0.6)**2+(y-0.35)**2)
    return np.maximum.reduce([bump1,bump2,bump3,np.zeros_like(bump1)])
def Volcano(x,y):
    r = np.sqrt((x-0.5)**2+(y-0.5)**2)
    volcano = 0.05+1.5*(1+x)*(r**2-6*r**4)
    return np.maximum(volcano, np.zeros_like(x))

def Volcano_gauss(x,y) :
    return(np.exp(-3*((x-0.5)**2+(y-0.5)**2))-0.7*np.exp(-9*((x-0.5)**2+(y-0.5)**2)))

def GenerateRHS(height,params):
    α,β,γ,h = params
    hx,hy = np.gradient(height,h)
    Intensity = (α*hx+β*hy+γ)/np.sqrt(1+hx**2+hy**2)
    Omega = height>0
    return Intensity,Omega

def Intensite_discontinue(height,masque,params):
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


def genere_intensite(H,params) :
    α,β,γ,h = params
    H2_x = np.roll(H,1,axis=1)
    H2_y = np.roll(H,1,axis=0)
    DXH = (H-H2_x)/h
    DYH = (H-H2_y)/h
    Intensity = (α*DXH+β*DYH+γ)/np.sqrt(1+DXH**2+DYH**2)
    return Intensity


# def gaussienne_x(x,y) :
#     return(np.exp(-17*(x-0.5)**2))
def gaussienne(x,y) :
    return(np.exp(-5*((x-0.5)**2+(y-0.5)**2)))

def gaussienne_x(x,y) :
    return(np.exp(-3*(1-x)**2))

def parametres(azimuth, elevation) :
    azimuth_rad = np.pi*azimuth/180
    elevation_rad = np.pi*elevation/180
    return(-np.cos(azimuth_rad)*np.cos(elevation_rad), -np.sin(azimuth_rad)*np.cos(elevation_rad),np.sin(elevation_rad))


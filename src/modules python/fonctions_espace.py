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

def GenerateRHS(height,params):
    α,β,γ,h = params
    hy,hx = np.gradient(height,h)
    Intensity = (α*hx+β*hy+γ)/np.sqrt(1+hx**2+hy**2)
    Omega = height>0
    return Intensity,Omega

def gaussienne_x(x,y) :
    return(np.exp(-10*(x-0.5)**2))

def gaussienne(x,y) :
    return(np.exp(-5*((x-0.5)**2+(y-0.5)**2)))

def parametres(azimuth, elevation) :
    azimuth_rad = np.pi*azimuth/180
    elevation_rad = np.pi*elevation/180
    return(-np.cos(azimuth_rad)*np.cos(elevation_rad), -np.sin(azimuth_rad)*np.cos(elevation_rad),np.sin(elevation_rad))

def GenerateIntDiscont(height,params) :
    (n,m) = height.shape
    α,β,γ,h = params
    hx,hy = np.gradient(height,h)
    hx[n//2,:] = np.zeros(m)
    hx[n//2-1,:] = np.zeros(m)
    Intensity = (α*hx+β*hy+γ)/np.sqrt(1+hx**2+hy**2)
    Omega = height>0
    return Intensity,Omega


def GenerateRHS_discontinue(height,masque,params):
    α,β,γ,h = params
    hx,hy = np.gradient(height,h)
    (n,m) = height.shape
    for i in range(n) :
        for j in range(m) :
            try :
                if masque[i,j] % 2 == 0 :
                    hx[i,j] = (height[i+1,j]-height[i,j])/h
                elif masque[i,j] % 3 == 0 :
                    if j == 0 :
                        hy[i,j] = 0
                    else :
                        hy[i,j] = (height[i,j]-height[i,j-1])/h
                elif masque[i,j] % 5 == 0 :
                    if i == 0 :
                        hx[i,j] == 0
                    else :
                        hx[i,j] = (height[i,j]-height[i-1,j])/h
                elif masque[i,j] % 7 == 0 :
                    hy[i,j] = (height[i,j+1]-height[i,j])/h
            except IndexError :
                if j == n-1 :
                    hy[i,j] = 0
                if i == m-1 :
                    hx[i,j] == 0

    Intensity = (α*hx+β*hy+γ)/np.sqrt(1+hx**2+hy**2)
    return Intensity



def GenerateIntIter(U,I,params) :
    (n,m) = U.shape
    α,β,γ,h = params
    p0 = np.array([-α,-β,γ])
    Intensity = np.zeros((n,m))
    hy,hx = np.gradient(U,h)

    for i in range(n) :
        for j in range(m) :
            p = p0*I[i,j]
            x = np.array([-hx[i,j],-hy[i,j],1])
            w= x
            k = x-p-np.dot(x-p,p)*p
            x = p + np.sqrt(1-np.dot(p,p))*((k)/np.sqrt(np.dot(k,k)))
            Intensity[i,j] = x[2]

    return(Intensity)

def IntInit(I,params) :
    n,m = I.shape
    I0 = np.zeros((n,m))
    alpha,beta,gamma,h = params
    p0 = np.array([-alpha,-beta,gamma])
    u = np.array([1,0,alpha/gamma])/np.sqrt(np.dot(np.array([1,0,alpha/gamma]),np.array([1,0,alpha/gamma])))
    v = np.array([-(beta*alpha)/(gamma**2+alpha**2),1,(beta/gamma) - (beta*alpha**2)/(gamma**3+alpha**2*gamma)])/np.sqrt(np.dot(np.array([-(beta*alpha)/(gamma**2+alpha**2),1,(beta/gamma) - (beta*alpha**2)/(gamma**3+alpha**2*gamma)]),np.array([-(beta*alpha)/(gamma**2+alpha**2),1,(beta/gamma) - (beta*alpha**2)/(gamma**3+alpha**2*gamma)])))
    for i in range(n) :
        for j in range(m) :
            a = 2*np.pi*np.random.rand()
            R = np.sqrt(1-I[i,j]**2)
            x = I[i,j]*p0 + R*(np.cos(a)*u + np.sin(a)*v)
            I0[i,j] = x[2]
    return(I0)

I = np.zeros((1,1))
params = 1,1,1,1

a = IntInit(I,params)








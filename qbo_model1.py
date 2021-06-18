#import os
#from netcdf4 import Dataset as netcdf_dataset
#from cartopy import config
#import cartopy.crs as ccrs
#import matplotlib.path as mpath
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# %% Parameters & Initializations

# beta plane approximation
omega = 7.29*10**(-5)               # rotation rate of earth (rad/s)
a = 6.371*10**6                     # radius of earth (m)
phi = 0                             # latitude for now
beta = 2*omega*np.cos(phi)*(1/a)    

# stepping & gridding
dt = 28800                          # time step (seconds); 1/3 day = 8 hrs
iters = 10000                       # number of time steps; 
t = np.arange(0,iters,dt)           # time grid

z0 = 17000                          # lower boundary (m); 17 km = tropical tropopause height (technically z(1/2))
zt = 70000                          # upper boundary (m); 70 km = mid-mesosphere 
z = np.linspace(z0,zt,121)          # vertical height grid; 120 levels (want these to be outer edges)
dz = np.zeros(len(z)-1)
for i in range(len(dz)):
    dz[i] = z[i+1] - z[i]           # layer heights

# other parameters for momentum equation
Kz = 0.3                            # vertical eddy diffusion coefficient (m^2/s)

# wave EP flux parameters
# Kelvin wave (F1)
F0k = 1                             # Plumb 77 flux init
ck = 25                             # phase speed in m/s  (Hamilton)
kk = 1                              # wavenumber  (Hamilton)

# Yanai (RG wave) (F2)
F0r = 1                             # Plumv 77 flux init
cr = -25                            # phase speed in m/s (negative for RG) (Hamilton)
kr = -4                             # wavenumber (neg for RG)  (Hamilton)

# Brunt Vaisala for EP flux WKB approximations
g = 9.81                            # gravity accel (m/s^2)
R = 287.058                         # specific gas constant (J/(kg*K)) 

# temperature and density profiles
T = pd.DataFrame.to_numpy(pd.read_csv('/Users/anncaseyhughes/Desktop/tropics_tempProfile.csv',engine='python',error_bad_lines=False))                   # import from Ananthayasam
dens = pd.DataFrame.to_numpy(pd.read_csv('/Users/anncaseyhughes/Desktop/tropics_densityProfile.csv',engine='python',error_bad_lines=False))

# create pressure profile
p = np.zeros(len(z))
sum = 0

for j in range(len(z)):
    if j == 0:
        sum = 0
    elif j == len(dz):
        sum += dz[j-1]/T[j]                     # deal with the fact that dz is midlayer defined and temperature is layer edge
    else:
        sum += dz[j]/T[j] 
    p[j] = 10000*np.exp((-g/R)*sum)             # 10000 Pa = 100 mb, height of tropical tropopause
    

# potential temperature profile
kappa = 0.286                       # exponent of pot temp equation (R/cp)
theta = np.zeros(len(T))            # potential temp profile
for j in range(len(T)):
    theta[j] = T[j]*(100000/(p[j]))**kappa                           # T(p0/p)^k
    
dThdz = np.zeros(len(z)-1)          # virtual potential temp lapse rate for BV freq (in layers)
for j in range(len(dThdz)):
    dThdz[j] = (theta[j+1]-theta[j])/(z[j+1]-z[j])
       
mTheta = np.zeros(len(dThdz))       # mean pot temp for midlayer BV freq computation
for j in range(len(dThdz)):
    mTheta[j] = (theta[j+1] + theta[j])/2
    
N = ((g/mTheta)*(dThdz))**0.5       # Brunt-Vaisala freq (theta version)


# Fels '81 Newtonian Cooling (just for CO2) INITIALIZATIONS (vertical wavenumber requires u and that needs to be calculated each time step)
theta0 = 960                         # basic state potential temperature for CO2  
zkm = z/1000                      
az = 0.422 + 0.001625*(zkm - 62.5) - 0.007125*(1 + (zkm - 62.5)**2)**0.5                 # a(z)
bz = 0.646 + 0.032*(zkm - 39.5) + 0.018*(9 + (zkm - 39.5)**2)**0.5                       # b(z)
cz = 0.045 + 0.000175*(zkm - 40)*((zkm - 40) - (1 + (zkm - 40)**2)**0.5)                 # c(z)

az = az[:-1]                        # trim so they're defined midlayer for cooling coefficient
bz = bz[:-1]
cz = cz[:-1]


# vertical wavenumber for wave momentum flux
#mk = np.zeros(len(z)-1)
#mr = np.zeros(len(z)-1)
#for j in range(len(mk)):
#    mk[j] = 0.000125                  # ave vertical wavenumber for Kelvin waves (Andrews 87)
#    mr[j] = 0.000167                  # "" Rossby-Gravity waves (Andrews 87)
#mk = 0.000125
#mr = 0.000167

# Newtonian cooling
#mT = np.zeros(len(z)-1)       # mean pot temp for midlayer BV freq computation
#for j in range(len(z)-1):
#    mT[j] = (T[j+1] + T[j])/2
#    mT[j] = (theta[j+1] + theta[j])/2

#alphak = np.zeros(len(z)-1)
#alphar = np.zeros(len(z)-1)
#for j in range(len(z)-1):
#    alphak[j] = (theta0/mT[j])**2 * np.exp((-theta0/mT[j]) * (az[j] + (bz[j]*mk[j]**2)/(cz[j]+mk[j]**1.5)))
#    alphar[j] = (theta0/mT[j])**2 * np.exp((-theta0/mT[j]) * (az[j] + (bz[j]*mr[j]**2)/(cz[j]+mr[j]**1.5)))
#    alphak[j] = (theta0/mT[j])**2 * np.exp((-theta0/mT[j]) * (az[j] + (bz[j]*mk**2)/(cz[j]+mk**1.5)))
#    alphar[j] = (theta0/mT[j])**2 * np.exp((-theta0/mT[j]) * (az[j] + (bz[j]*mr**2)/(cz[j]+mr**1.5)))
alpha = 1.5 + (10**(-6))*np.tanh((zkm[:-1]-35)/7)                            # from Holton & Mass 1976

# %% Numerical integration of vertical EP flux 

# initialize u vector at all heights
u = np.zeros(len(z)-1)              # init zonal mean zonal wind
diagwind = []
diagpts = [1,5,10]

#for i in range(iters):
for i in range(2):
#    if i in diagpts:
#        diagwind = np.append(diagwind,u)
#    ubar = np.mean(u)    
#    
#    # vertical wavenumber for wave momentum flux generated by the zero frequency nuik
##    mk = np.abs(N/(ck-ubar))                            # dispersion relation from Ryu 2007
##    mr = np.abs((N/(ck-ubar))+((N*beta)/((kr^2)*(cr-ubar))))                                       # dispersion relation from Appendix A of Smyth 2007
#    mk = np.zeros(len(z)-1)
#    mr = np.zeros(len(z)-1)
#    for j in range(len(mk)):
#        mk[j] = 0.000125                  # ave vertical wavenumber for Kelvin waves (Andrews 87)
#        mr[j] = 0.000167                  # "" Rossby-Gravity waves (Andrews 87)
#
#    # Newtonian cooling
#    mT = np.zeros(len(z)-1)       # mean pot temp for midlayer BV freq computation
#    for j in range(len(z)-1):
#        mT[j] = (T[j+1] + T[j])/2
#    
#    alphak = np.zeros(len(z)-1)
#    alphar = np.zeros(len(z)-1)
#    for j in range(len(z)-1):
#        alphak[j] = (theta0/mT[j])**2 * np.exp((-theta0/mT[j]) * (az[j] + (bz[j]*mk[j]**2)/(cz[j]+mk[j]**1.5)))
#        alphar[j] = (theta0/mT[j])**2 * np.exp((-theta0/mT[j]) * (az[j] + (bz[j]*mr[j]**2)/(cz[j]+mr[j]**1.5)))
     
    # Flux functions for wave momentum
#    gk = (alphak*N)/(kk*(ck-ubar)**2)                                               # inner function for flux of Kelvin wave momentum
#    gr = ((N*beta*alphar)/(kr**3*(cr-ubar)**3))*((1+(cr-ubar))/(beta*kr**2))           # inner function for flux of mixed RG wave momentum
#    gk = (alpha*N)/(kk*(ck-ubar)**2)                                               # inner function for flux of Kelvin wave momentum
#    gr = ((N*beta*alpha)/(kr**3*(cr-ubar)**3))*((1+(cr-ubar))/(beta*kr**2))           # inner function for flux of mixed RG wave momentum
    gk = np.zeros(len(z)-1)
    gr = np.zeros(len(z)-1)
    for j in range(len(gk)):
        gk[j] = (alpha[j]*N[j])/(kk*(ck-u[j])**2)
        gr[j] = ((N[j]*beta*alpha[j])/(kr**3*(cr-u[j])**3))*((1+cr-u[j])/(beta*kr**2))
        
    Fk = np.zeros(len(z))               # initialize momentum flux for Kelvin waves
    Fr = np.zeros(len(z))               # initialize momentum flux for mixed RG 
    Fk[0] = 1
    Fr[0] =1
    for j in range(len(z)):
        sum1 = 0
        sum2 = 0
        for k in range(j):
            sum1 += dz[k]*gk[k]
            sum2 += dz[k]*gr[k]
        Fk[j] = Fk[0] * np.exp(-sum1)
        Fr[j] = Fr[0] * np.exp(-sum2)
        
    
    # forward difference flux vectors to get convergence terms
    kelv_fluxconv = np.zeros(len(z)-1)
    mrg_fluxconv = np.zeros(len(z)-1)
    for j in range(len(kelv_fluxconv)):
        kelv_fluxconv[j] = (-1/dens[j])*(Fk[j+1]-Fk[j])/dz[j]
        mrg_fluxconv[j] = (-1/dens[j])*(Fr[j+1]-Fr[j])/dz[j]
        
    # create eddy diffusion vectors
    eddydiff = np.zeros(len(z)-1)
    for j in range(len(z)-3):
        eddydiff[j+1] = Kz*(u[j+1]-2*u[j]-u[j])/(dz[j-1]**2)
    eddydiff[0] = eddydiff[1]
    eddydiff[len(z)-2] = eddydiff[len(z)-3]    
    
    # create temporal tendency update vector    
    h = kelv_fluxconv + eddydiff + mrg_fluxconv
    hk = kelv_fluxconv + eddydiff
    hr = mrg_fluxconv + eddydiff
    
    for j in range(len(u)):
        u[j] = u[j] +dt*h[j]

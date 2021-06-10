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
kk = 1                             # wavenumber  (Hamilton)

# Yanai (RG wave) (F2)
F0r = 1                             # Plumv 77 flux init
cr = -25                            # phase speed in m/s (negative for RG) (Hamilton)
kr = -4                             # wavenumber (neg for RG)  (Hamilton)

# Brunt Vaisala for EP flux WKB approximations
g = 9.81                            # gravity accel (m/s^2)

# temperature and density profiles
T = pd.DataFrame.to_numpy(pd.read_csv('/Users/anncaseyhughes/Desktop/tropics_tempProfile.csv',engine='python',error_bad_lines=False))                   # import from Ananthayasam
dens = pd.DataFrame.to_numpy(pd.read_csv('/Users/anncaseyhughes/Desktop/tropics_densityProfile.csv',engine='python',error_bad_lines=False))

# potential temperature profile
kappa = 0.286                       # exponent of pot temp equation (R/cp)
R = 287.058                         # specific gas constant (J/(kg*K)) 
theta = np.zeros(len(T))            # potential temp profile
for i in range(len(T)):
    theta[i] = T[i]*(100000/(dens[i]*R*T[i]))**kappa                           # T(p0/p)^k
    
thetaV = np.mean(theta)             # ambient potential temperature for BV freq (on layer EDGES)
dThdz = np.zeros(len(z)-1)          # virtual potential temp lapse rate for BV freq (in layers)
for i in range(len(dThdz)):
    dThdz[i] = (theta[i+1]-theta[i])/(z[i+1]-z[i])
    
zero_dThdz = np.zeros(len(dThdz))
for i in range(len(dThdz)):
    if dThdz[i] < 0:
        zero_dThdz[i] = 0
    
N = ((g/thetaV)*(zero_dThdz))**0.5       # Brunt-Vaisala freq (theta version)

# Fels '81 Newtonian Cooling (just for CO2)
theta0 = 960                         # basic state potential temperature for CO2
T0 = np.mean(T)                      # basic state temperature  
zkm = z/1000                      
az = 0.422 + 0.001625*(zkm - 62.5) - 0.007125*(1 + (zkm - 62.5)**2)**0.5                 # a(z)
bz = 0.646 + 0.032*(zkm - 39.5) + 0.018*(9 + (zkm - 39.5)**2)**0.5                       # b(z)
cz = 0.045 + 0.000175*(zkm - 40)*((zkm - 40) - (1 + (zkm - 40)**2)**0.5)                 # c(z)

#alpha = (theta0/T)**2 * np.exp(-theta0/T) * (az + (bz*m**2)/())

# %% Numerical integration of vertical EP flux 

# initialize u vector at all heights

u = np.zeros(len(z)-1)              # init zonal mean zonal wind
#for i in range(len(t)):
    # create flux vectors
    
m = np.zeros(len(z))
alpha = np.zeros((len(z)-1,len(m)))
for j in range(len(z)-1):
    if u[j] == 0:
        m[j] = N[j]/1
    else: 
        m[j] = N[j]/u[j]                  # N/U (Beres 2004)

m[len(z)-1] = m[len(z)-2]

for i in range(len(z)-1):
    alpha
alpha = (theta0/T)**2 * np.exp(-theta0/T) * (az + (bz*m**2)/(cz+m**1.5))

g1 = (alpha[:,0]*N)/(kk*(ck-u)**2)                                               # inner function for flux of Kelvin wave momentum
g2 = ((N*beta*alpha[:,0])/(kr**3*(cr-u)**3))*((1+(cr-u))/(beta*kr**2))           # inner function for flux of mixed RG wave momentum
# [:,0] above for just wavenumber 1

F1 = np.zeros(len(z))               # initialize momentum flux for Kelvin waves
F2 = np.zeros(len(z))               # initialize momentum flux for mixed RG waves
for j in range(len(z)):
    sum1 = 0
    sum2 = 0
    for k in range(j):
        sum1 += dz[k]*g1[k]
        sum2 += dz[k]*g2[k]
    F1[j] = F0k * np.exp(-sum1)
    F2[j] = F0r * np.exp(-sum2)

# forward difference flux vectors to get convergence terms
kelv_fluxconv = np.zeros(len(z)-1)
mrg_fluxconv = np.zeros(len(z)-1)
for j in range(len(kelv_fluxconv)):
    kelv_fluxconv[j] = (-1/dens[j])*(F1[j+1]-F1[j])/dz[j]
    mrg_fluxconv[j] = (-1/dens[j])*(F2[j+1]-F2[j])/dz[j]
    
# create eddy diffusion vectors
eddydiff = np.zeros(len(z)-1)
for j in range(len(z)-3):
    eddydiff[j+1] = Kz*(u[j+1]-2*u[j]-u[j])/(dz[j-1]**2)
eddydiff[0] = eddydiff[1]
eddydiff[len(z)-2] = eddydiff[len(z)-3]    
    
h = kelv_fluxconv + eddydiff + mrg_fluxconv
for j in range(len(u)):
    u[j] = u[j] +dt*h[j]
    






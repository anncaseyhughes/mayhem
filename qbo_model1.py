#import os
#from netcdf4 import Dataset as netcdf_dataset
#from cartopy import config
#import cartopy.crs as ccrs
#import matplotlib.path as mpath
#import matplotlib.pyplot as plt
import numpy as np

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

z0 = 17000                          # lower boundary (m); 17 km = tropical tropopause height
zt = 70000                          # upper boundary (m); 70 km = mid-mesosphere
z = np.linspace(z0,zt,120)          # vertical height grid; 120 levels

# other parameters for momentum equation
Kz = 0.3                            # vertical eddy diffusion coefficient (m^2/s)
# rho0 = ??                           # basic state density (kg/m^3)

# wave EP flux parameters
# Kelvin wave (F1)
# F0k = ???
ck = 25                             # phase speed in m/s
kk = -1                             # wavenumber

# Yanai (RG wave) (F2)
# F0r = ???
cr = -25                            # phase speed in m/s (negative for RG)
kr = -4                             # wavenumber (neg for RG)

# Fels '81 Newtonian Cooling (just for CO2)
theta0 = 960
#T0 = ???                            # basic state temperature?                        
# beta(z) = a(z) = 0.422 + 0.001625*(z - 62.5) - 0.007125*(1 + (z - 62.5)**2)**0.5
# gamma(z) = b(z) = 0.646 + 0.032*(z - 39.5) + 0.018*(9 - (z - 39.5)**2)**0.5
# delta(z) = c(z) = 0.045 + 0.000175*(z - 40)*((z - 40) - (1 + (z - 40)**2)**0.5)
# alpha(z,m) = 

# Brunt Vaisala for EP flux WKB approximations
g = 9.8                             # gravity accel (m/s^2)
#thetaV = ?????                       # ambient vritual potential temperature for BV freq
#dThdz = ????                         # virtual potential temp lapse rate for BV freq
#N = ((g/theta0)*(dThdz))**0.5              # Brunt-Vaisala freq

# %% Numerical integration of vertical EP flux 

# initialize flux vectors
Fk = np.zeros(len())
#Fk = F0k * exp(PLACE TRAPEZOID RULE QUADRATURE HERE)
#Fr = F0r + exp(PLACE TRAPEZOID RULE QUADRATURE HERE)





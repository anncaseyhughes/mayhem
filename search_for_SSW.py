# import os
from netCDF4 import Dataset
# from cartopy import config
import cartopy.crs as ccrs
import matplotlib.path as mpath
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

#dataset = xr.open_dataset('/Volumes/pond/mark/subset-B_2000_CAM5_nlev53-0001.nc')
#dataset = xr.open_dataset('/Volumes/pond/mark/subset-B_2000_CAM5_nlev53-0002.nc')
#dataset = xr.open_dataset('/Volumes/pond/mark/subset-B_2000_CAM5_nlev53-0003.nc')
#dataset = xr.open_dataset('/Volumes/pond/mark/subset-B_2000_CAM5_nlev53-0004.nc')
dataset = xr.open_dataset('/Volumes/pond/mark/subset-B_2000_CAM5_nlev53-0005.nc')

T = dataset.variables['T'][:, 37, :, :]			# just want north pole

#go through all time, all space, 




# Plotting
# lats = dataset.variables['lat'][:]
# lons = dataset.variables['lon'][:]

## theta = np.linspace(0, 2*np.pi, 100)
# center, radius = [0.5, 0.5], 0.5
# verts = np.vstack([np.sin(theta), np.cos(theta)]).T
# circle = mpath.Path(verts * radius + center)
# #
# fig = plt.figure()
# ax = plt.axes(projection=ccrs.NorthPolarStereo())
# ##plt.title('North Polar Stereographic Projection of Wind Speed at 100 hPa Isobar')
# ##plt.title('North Polar Stereographic Projection of Wind Speed at 0 hPa Isobar')
# plt.title('North Polar Stereographic Projection of Temperature at 100 hPa Isobar (dataset 5)')
# ##plt.title('North Polar Stereographic Projection of Temperature at 0 hPa Isobar')
# #
# #
# ax.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())
# ax.set_boundary(circle, transform=ax.transAxes)
# 
# ##plt.contourf(lons, lats, U, 60, transform=ccrs.PlateCarree())
# ##filled_c = plt.contourf(lons, lats, U, 60, transform=ccrs.PlateCarree())
# plt.contourf(lons, lats, T, 60, transform=ccrs.PlateCarree())
# filled_c = plt.contourf(lons, lats, T, 60, transform=ccrs.PlateCarree())
# #
# ax.coastlines()
# ax.gridlines()
# #
# cbar = fig.colorbar(filled_c, orientation='horizontal')
# cbar.set_label('degrees in K')
# #
# ##plt.savefig('northpolarstereo.eps', format='eps')
# ##plt.savefig('temp_0hPa.eps', format='eps')
# #plt.savefig('temp_100hPa.eps', format='eps')
# plt.savefig('newset_temptest5.eps',format='eps')
# plt.show()
#

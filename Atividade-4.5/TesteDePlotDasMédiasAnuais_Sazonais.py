import xarray as xr
import numpy as np
import geopandas as gpd
from PIL import Image

import matplotlib.pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
import matplotlib.gridspec as gridspec

# Leitura dos estados do shapefile
estados = gpd.read_file('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/BR_UF_2022.shp')

# Leitura do conjunto de dados de precipitação
dset = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/Chirps1AmericaDoSul1991_2020.nc')

dsetMediaAnual = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/Chirps1AmericaDoSul1991_2020MediaAnual.nc')

dsetMediaSazonal = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/Chirps1AmericaDoSul1991_2020MediaSazonal.nc')

# print(dset['time.season'].values)
# print(dset['time.month'].values)


lat = dset.latitude.values
lon = dset.longitude.values
var = dset['precip']

# latMediaAnual = dsetMediaAnual.latitude.values
# lonMediaAnual = dsetMediaAnual.longitude.values
# varMediaAnual = dsetMediaAnual['precip']

# latMediaSazonal = dsetMediaSazonal.latitude.values
# lonMediaSazonal = dsetMediaSazonal.longitude.values
# varMediaSazonal = dsetMediaSazonal['precip']

# # Calcula a média anual dos dados de precipitação
clim = np.mean(var, axis=0)
# climMediaAnual = np.mean(varMediaAnual, axis=0)
# climMediaSazonal = np.mean(varMediaSazonal, axis=0)
# print(clim.coords)
# print(climMediaAnual.coords)
# print(climMediaSazonal.coords)

sclim = var.groupby('time.season').mean('time')
# sclimMediaAnual = varMediaAnual.groupby('time.season').mean('time')
# sclimMediaSazonal = varMediaSazonal.groupby('time.season').mean('time')

# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.set_extent((-90.0, -30.0, -60.0, 15.0))
# ax.coastlines(resolution='110m', color = 'black')
# ax.add_feature(cfeature.BORDERS, linewith= 0.5, edgecolor = 'black')

# # Configuração das grades do mapa
# gl = ax.gridlines(draw_labels=True)
# gl.xlines = False
# gl.ylines = False
# gl.right_labels = False
# gl.top_labels = False

# levels = np.linspace(0.0, 30.0,11)

def create_plot(ax,nrow, ncol, data, tlon, season = ''):
    ax.set_extent((-90.0, -30.0, -60.0, 15.0))
    ax.coastlines(resolution='110m', color = 'black')
    ax.add_feature(cfeature.BORDERS, linewith= 0.5, edgecolor = 'black')

    gl = ax.gridlines(draw_labels=True)
    gl.xlines = False
    gl.ylines = False
    gl.right_labels = False
    gl.top_labels = False

    levels = np.linspace(0.0, 18.0, 11)

    cnplot = ax.contourf(data.longitude, data.latitude, data.sel(season=season).values, levels=levels, cmap='YlGnBu', extend='max')
    ax.text(tlon, -55, season, bbox=dict(facecolor='white', alpha=0.7))

    return cnplot

# def create_plot_MediaAnual(ax,nrow, ncol, data, tlon, season = ''):
#     ax.set_extent((-90.0, -30.0, -60.0, 15.0))
#     ax.coastlines(resolution='110m', color = 'black')
#     ax.add_feature(cfeature.BORDERS, linewith= 0.5, edgecolor = 'black')

#     gl = ax.gridlines(draw_labels=True)
#     gl.xlines = False
#     gl.ylines = False
#     gl.right_labels = False
#     gl.top_labels = False

#     levels = np.linspace(0.0, 18.0, 11)

#     cnplotMediaAnual = ax.contourf(data.longitude, data.latitude, data.sel(season=season).values, levels=levels, cmap='YlGnBu', extend='max')
#     ax.text(tlon, -55, season, bbox=dict(facecolor='white', alpha=0.7))

#     return cnplotMediaAnual

# def create_plot_MediaSazonal(ax,nrow, ncol, data, tlon, season = ''):
#     ax.set_extent((-90.0, -30.0, -60.0, 15.0))
#     ax.coastlines(resolution='110m', color = 'black')
#     ax.add_feature(cfeature.BORDERS, linewith= 0.5, edgecolor = 'black')

#     gl = ax.gridlines(draw_labels=True)
#     gl.xlines = False
#     gl.ylines = False
#     gl.right_labels = False
#     gl.top_labels = False

#     levels = np.linspace(0.0, 18.0, 11)

#     cnplotMediaSazonal = ax.contourf(data.longitude, data.latitude, data.sel(season=season).values, levels=levels, cmap='YlGnBu', extend='max')
#     ax.text(tlon, -55, season, bbox=dict(facecolor='white', alpha=0.7))

#     return cnplotMediaSazonal

fig = plt.figure()

gs = gridspec.GridSpec(4, 3, hspace=0.12, wspace=0.00)

# Infer seasons from time values
# months = dset['time.month'].values
# seasons = ['DJF', 'MAM', 'JJA', 'SON']
# season_labels = [seasons[(month-1)//3] for month in months]


seasons = ['DJF', 'MAM', 'JJA', 'SON']

# for i, season in enumerate(season_labels):
#     row = i // 4  # Use integer division to get the row index
#     col = i % 4  # Use modulo to get the column index
create_plot(plt.subplot( projection=ccrs.PlateCarree()), sclim, -44, seasons)


    # create_plot_MediaAnual(plt.subplot(gs[row, col], projection=ccrs.PlateCarree()), row, col, sclimMediaAnual, -44, season)
    # create_plot_MediaSazonal(plt.subplot(gs[row, col], projection=ccrs.PlateCarree()), row, col, sclimMediaSazonal, -44, season)

cax = plt.axes([0.2, 0.065, 0.6, 0.069])
cbar = plt.colorbar(create_plot(plt.subplot(gs[0, 0], projection=ccrs.PlateCarree()), 0, 0, sclim, -44, 'DJF'), cax=cax, orientation='horizontal', pad=0.4)
# cbar = plt.colorbar(create_plot_MediaAnual(plt.subplot(gs[0, 1], projection=ccrs.PlateCarree()), 0, 0, sclimMediaAnual, -44, 'DJF'), cax=cax, orientation='horizontal', pad=0.4)
# cbar = plt.colorbar(create_plot_MediaSazonal(plt.subplot(gs[0, 2], projection=ccrs.PlateCarree()), 0, 0, sclimMediaSazonal, -44, 'DJF'), cax=cax, orientation='horizontal', pad=0.4)

plt.show()
import xarray as xr
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
import matplotlib.gridspec as gridspec

# Leitura dos estados do shapefile
estados = gpd.read_file('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/BR_UF_2022.shp')

# Leitura do conjunto de dados de precipitação
dsetMediaAnual = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/Chirps1AmericaDoSul1991_2020MediaAnual.nc')

# Extrai as coordenadas de latitude e longitude do arquivo de dados
lat = dsetMediaAnual.latitude.values
lon = dsetMediaAnual.longitude.values

# Seleciona a variável 'precip' do arquivo de dados
var = dsetMediaAnual['precip']

# Calcula a média sazonal dos dados de precipitação
sclim = var.groupby('time.season').mean('time')

# Obtém as primeiras 9 temporadas 'JJA'
primeiras_jja = dsetMediaAnual['time.season'][:9]

# 3 - Plots
def create_plot(ax, data, tlon, season=''):
    ax.set_extent((-90.0, -30.0, -60.0, 15.0))
    ax.coastlines(resolution='110m', color='black')
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')

    gl = ax.gridlines(draw_labels=True)
    gl.xlines = False
    gl.ylines = False
    gl.right_labels = False
    gl.top_labels = False

    levels = np.linspace(0.0, 400.0, 11)
    cnplot = ax.contourf(data.longitude, data.latitude, data.sel(season=season), levels=levels, cmap='YlGnBu', extend='max')
    ax.text(tlon, -55, season, bbox=dict(facecolor='white', alpha=0.7))

    return cnplot

# Cria uma figura e uma grade de subplots
fig = plt.figure(figsize=(15, 15))
gs = gridspec.GridSpec(3, 3, hspace=0.4, wspace=0.4)

# Cria e plota os subplots para cada temporada
for i, season in enumerate(primeiras_jja):
    row = i // 3
    col = i % 3
    create_plot(plt.subplot(gs[row, col], projection=ccrs.PlateCarree()), sclim, -44, season)

# Adiciona uma barra de cores
cax = plt.axes([0.2, 0.065, 0.6, 0.02])
cbar = plt.colorbar(create_plot(plt.subplot(gs[0, 0], projection=ccrs.PlateCarree()), sclim, -44, primeiras_jja[0]), cax=cax, orientation='horizontal', pad=0.4)

# Salva a figura e a exibe
plt.savefig('s2-CHIRPSMediaAnual.png', dpi=300)
plt.show()

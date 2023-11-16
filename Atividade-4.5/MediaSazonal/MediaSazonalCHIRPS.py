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
dsetMediaSazonal = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/Chirps1AmericaDoSul1991_2020MediaSazonal.nc')
# print(dset)
# print(dset['time_bnds'])


# # Extrai as coordenadas de latitude e longitude do arquivo de dados
lat = dsetMediaSazonal.latitude.values
lon = dsetMediaSazonal.longitude.values

# Seleciona a variável 'precip' do arquivo de dados
var = dsetMediaSazonal['precip']

# # Calcula a média anual dos dados de precipitação
clim = np.mean(var, axis=0)
# # print(clim.coords)

# # Calcula a média sazonal dos dados de precipitação
sclim = var.groupby('time.season').mean('time')

# # print(sclim)
# # print(sclim['season'].values)

# #-------------------------------------------------------------------------#
# # 3 - Plots

# Configuração da projeção do mapa e sua extensão geográfica
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent((-90.0, -30.0, -60.0, 15.0))
ax.coastlines(resolution='110m', color='black')
ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')

# Configuração das grades do mapa
gl = ax.gridlines(draw_labels=True)
gl.xlines = False
gl.ylines = False
gl.right_labels = False
gl.top_labels = False

# Define os níveis de contorno
levels = np.linspace(0.0, 500.0, 11)

# Inicializa um gráfico de contorno preenchido com os dados de precipitação (clim)
cnplot = ax.contourf(clim.longitude, clim.latitude, clim, cmap='YlGnBu',levels = levels, extend='max')

# Adiciona uma barra de cores (colorbar) na parte inferior do gráfico
cbar = plt.colorbar(cnplot, orientation='horizontal', pad=0.07, shrink=0.6)
cbar.set_label('precipitation (mm/day) \n CHIRPS(2002)')

# Define o título do gráfico
ax.set_title('(matplotlib)')

# Plota os limites dos estados (supondo que sejam limites de estados geográficos)
estados.plot(ax=ax, color='none', edgecolor='black')

# Salva a figura em um arquivo chamado 's1-matplot.png' com uma resolução de 300 dpi
plt.savefig('CHIRPSMediaSazonal.png', dpi=300)

# Exibe o gráfico na tela
plt.show()

def create_plot(ax, nrow, ncol, data, tlon, season=''):
    ax.set_extent((-90.0, -30.0, -60.0, 15.0))  # Define a extensão geográfica
    ax.coastlines(resolution='110m', color='black')  # Adiciona a linha da costa
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')  # Adiciona as fronteiras

    gl = ax.gridlines(draw_labels=True)
    gl.xlines = False
    gl.ylines = False
    gl.right_labels = False
    gl.top_labels = False

    levels = np.linspace(0.0, 500.0, 11)  # Define os níveis de contorno Aqui, você está criando uma lista de 11 valores espaçados uniformemente entre 0 e 100. Esses valores representam os níveis de contorno que serão usados no gráfico.
    cnplot = ax.contourf(data.longitude, data.latitude, data.sel(season=season), levels = levels, cmap='YlGnBu', extend='max')  # Cria um plot de contorno
    ax.text(tlon, -55, season, bbox=dict(facecolor='white', alpha=0.7))  # Adiciona o nome da estação

    return cnplot

# Cria uma figura e uma grade de subplots
fig = plt.figure()
# plt.show()
gs = gridspec.GridSpec(2, 2, hspace=0.12, wspace=0.00)

# Define as estações que deseja plotar
seasons = ['DJF', 'MAM', 'JJA', 'SON']

# Cria e plota os subplots para cada estação
for i, season in enumerate(seasons):
    row = i // 2
    col = i % 2
    create_plot(plt.subplot(gs[row, col], projection=ccrs.PlateCarree()), row, col, sclim, -44, season)

# Adiciona uma barra de cores
cax = plt.axes([0.2, 0.065, 0.6, 0.069])
cbar = plt.colorbar(create_plot(plt.subplot(gs[0, 0], projection=ccrs.PlateCarree()), 0, 0, sclim, -44, 'DJF'), cax=cax, orientation='horizontal', pad=0.4)

# Salva a figura e a exibe
plt.savefig('s2-CHIRPSMediaSazonal.png', dpi=300)
plt.show()
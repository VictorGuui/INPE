
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

latMediaAnual = dsetMediaAnual.latitude.values
lonMediaAnual = dsetMediaAnual.longitude.values
varMediaAnual = dsetMediaAnual['precip']

latMediaSazonal = dsetMediaSazonal.latitude.values
lonMediaSazonal = dsetMediaSazonal.longitude.values
varMediaSazonal = dsetMediaSazonal['precip']

# # Calcula a média anual dos dados de precipitação
clim = np.mean(var, axis=0)
climMediaAnual = np.mean(varMediaAnual, axis=0)
climMediaSazonal = np.mean(varMediaSazonal, axis=0)
# print(clim.coords)
# print(climMediaAnual.coords)
# print(climMediaSazonal.coords)

sclim = var.groupby('time.season').mean('time')
sclimMediaAnual = varMediaAnual.groupby('time.season').mean('time')
sclimMediaSazonal = varMediaSazonal.groupby('time.season').mean('time')




# # 3 - Plots
# print('Chirps')

# # Configuração da projeção do mapa e sua extensão geográfica
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.set_extent((-90.0, -30.0, -60.0, 15.0))
# ax.coastlines(resolution='110m', color='black')
# ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')

# # Configuração das grades do mapa
# gl = ax.gridlines(draw_labels=True)
# gl.xlines = False
# gl.ylines = False
# gl.right_labels = False
# gl.top_labels = False

# # Define os níveis de contorno
# levels = np.linspace(0.0, 500.0, 11)

# # Inicializa um gráfico de contorno preenchido com os dados de precipitação (clim)
# cnplot = ax.contourf(clim.longitude, clim.latitude, clim, cmap='YlGnBu', extend='max')

# # Adiciona uma barra de cores (colorbar) na parte inferior do gráfico
# cbar = plt.colorbar(cnplot, orientation='horizontal', pad=0.07, shrink=0.6)
# cbar.set_label('precipitation (mm/day) \n CHIRPS')

# # Define o título do gráfico
# ax.set_title(' (matplotlib)')

# # Plota os limites dos estados (supondo que sejam limites de estados geográficos)
# estados.plot(ax=ax, color='none', edgecolor='black')

# # Salva a figura em um arquivo chamado 's1-matplot.png' com uma resolução de 300 dpi
# plt.savefig('s1-CHIRPS.png', dpi=300)

# # Exibe o gráfico na tela
# plt.show()

print('Media Anual')

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
cnplotAnual = ax.contourf(climMediaAnual.longitude, climMediaAnual.latitude, climMediaAnual, cmap='YlGnBu', extend='max')

# Adiciona uma barra de cores (colorbar) na parte inferior do gráfico
cbarAnual = plt.colorbar(cnplotAnual, orientation='horizontal', pad=0.07, shrink=0.6)
cbarAnual.set_label('precipitation (mm/day) \n CHIRPS_Anual')

# Define o título do gráfico
ax.set_title(' (matplotlib)')

# Plota os limites dos estados (supondo que sejam limites de estados geográficos)
estados.plot(ax=ax, color='none', edgecolor='black')

# Salva a figura em um arquivo chamado 's1-matplot.png' com uma resolução de 300 dpi
plt.savefig('s1-CHIRPS_Anual.png', dpi=300)

# Exibe o gráfico na tela
plt.show()

# print('Media Sazonal')
# # Configuração da projeção do mapa e sua extensão geográfica
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.set_extent((-90.0, -30.0, -60.0, 15.0))
# ax.coastlines(resolution='110m', color='black')
# ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')

# # Configuração das grades do mapa
# gl = ax.gridlines(draw_labels=True)
# gl.xlines = False
# gl.ylines = False
# gl.right_labels = False
# gl.top_labels = False

# # Define os níveis de contorno
# levels = np.linspace(0.0, 500.0, 11)

# # Inicializa um gráfico de contorno preenchido com os dados de precipitação (clim)
# cnplotSazonal = ax.contourf(climMediaSazonal.longitude, climMediaSazonal.latitude, climMediaSazonal, cmap='YlGnBu', extend='max')

# # Adiciona uma barra de cores (colorbar) na parte inferior do gráfico
# cbarSazonal = plt.colorbar(cnplotSazonal, orientation='horizontal', pad=0.07, shrink=0.6)
# cbarSazonal.set_label('precipitation (mm/day) \n CHIRPSSazonal')

# # Define o título do gráfico
# ax.set_title(' (matplotlib)')

# # Plota os limites dos estados (supondo que sejam limites de estados geográficos)
# estados.plot(ax=ax, color='none', edgecolor='black')

# # Salva a figura em um arquivo chamado 's1-matplot.png' com uma resolução de 300 dpi
# plt.savefig('s1-CHIRPSSazonal.png', dpi=300)

# # Exibe o gráfico na tela
# plt.show()

# Abre as imagens salvas
img = Image.open("s1-CHIRPS.png")
img1 = Image.open("s1-CHIRPS_Anual.png")
img2 = Image.open("s1-CHIRPSSazonal.png")

# Redimensiona as imagens
img_size = img.resize((500, 768))
img1_size = img1.resize((500, 768))
img2_size = img2.resize((500, 768))

# Cria uma nova imagem branca
img3 = Image.new("RGB", (1500, 800), "white")

# Cola as imagens na nova imagem
img3.paste(img_size, (0, 0))
img3.paste(img1_size, (500, 0))
img3.paste(img2_size, (1000, 0))

# Salva a nova imagem
img3.save('s1-chirps-medias-anual-sazonal.png')
    

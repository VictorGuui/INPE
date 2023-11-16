import xarray as xr
import numpy as np
import geopandas as gpd
from PIL import Image
import matplotlib.pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
import matplotlib.gridspec as gridspec

# Carrega os dados geoespaciais dos estados do Brasil
estados = gpd.read_file('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/BR_UF_2022.shp')

# Carrega os dados climatológicos de precipitação
dset = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/arquivoClimatologia_1991_2020.nc')
lat = dset.lat.values
lon = dset.lon.values
var = dset['precip']

# Calcula a média anual dos dados de precipitação
clim = np.mean(var, axis=0)

# Calcula a média sazonal dos dados de precipitação
sclim = var.groupby('time.season').mean('time')

# Define uma função para criar subplots no mapa
def create_plot(ax, nrow, ncol, data, tlon, season=''):
    ax.set_extent((-90.0, -30.0, -60.0, 15.0))
    ax.coastlines(resolution='110m', color='black')
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')

    gl = ax.gridlines(draw_labels=True)
    gl.xlines = False
    gl.ylines = False
    gl.right_labels = False
    gl.top_labels = False

    levels = np.linspace(0.0, 10.0, 11)
    cnplot = ax.contourf(data.lon, data.lat, data.sel(season=season), cmap='YlGnBu', extend='max')
    ax.text(tlon, -55, season, bbox=dict(facecolor='white', alpha=0.7))

    return cnplot

# Cria uma figura e uma grade de subplots
fig = plt.figure()
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

# Salva a figura
plt.savefig('s2-matplot.png', dpi=300)
plt.show()

# Abre as imagens salvas
img = Image.open("s1-matplot.png")
img1 = Image.open("s2-matplot.png")

# Redimensiona as imagens
img_size = img.resize((670, 1024))
img1_size = img1.resize((670, 1024))

# Cria uma nova imagem branca
img2 = Image.new("RGB", (1340, 1024), "white")

# Cola as imagens na nova imagem
img2.paste(img_size, (0, 0))
img2.paste(img1_size, (670, 0))

# Salva a nova imagem
img2.save('s1-medias-anual-sazonal.png')

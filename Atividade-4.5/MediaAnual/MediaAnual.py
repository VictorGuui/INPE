import xarray as xr
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
import matplotlib.gridspec as gridspec
from PIL import Image

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
primeiras_jja = dsetMediaAnual['time.season'][:9]

# Seleciona os dados de precipitação para os primeiros 9 anos completos
dados_primeiros_anos = dsetMediaAnual['precip'].sel(time=slice('1991-06-01', '1999-08-31'))

# Calcula a média sazonal dos dados de precipitação para os primeiros 9 anos completos
sclim = var.groupby('time.season').mean('time')

clim = np.mean(var, axis=0)

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
levels = np.linspace(0.0, 550.0, 11)

# Inicializa um gráfico de contorno preenchido com os dados de precipitação (clim)
cnplot = ax.contourf(clim.longitude, clim.latitude, clim,levels=levels, cmap='YlGnBu', extend='max')

# Adiciona uma barra de cores (colorbar) na parte inferior do gráfico
cbar = plt.colorbar(cnplot, orientation='horizontal', pad=0.07, shrink=0.6)
cbar.set_label('precipitation (mm/day) \n CHIRPS Média Anual')

# Plota os limites dos estados (supondo que sejam limites de estados geográficos)
estados.plot(ax=ax, color='none', edgecolor='black')

# Salva a figura em um arquivo chamado 's1-matplot.png' com uma resolução de 300 dpi
plt.savefig('Media dos anos(CHIRPS Media Anual).png', dpi=300)

# Exibe o gráfico na tela
plt.show()

# 3 - Plots
def create_plot(ax, data, tlon, season, year=''):

    # A função é definida com o nome create_plot e aceita vários parâmetros, incluindo ax (o objeto de eixo no qual o gráfico será desenhado), data (os dados a serem plotados), tlon (um valor para a longitude), season (a estação do ano) e year (o ano para adicionar ao título).

    ax.set_extent((-90.0, -30.0, -60.0, 15.0))
    ax.coastlines(resolution='110m', color='black')
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')

    # Essas linhas configuram o mapa no objeto de eixo (ax). Elas definem a extensão geográfica, adicionam linhas costeiras e fronteiras ao mapa.

    gl = ax.gridlines(draw_labels=True)
    gl.xlines = False
    gl.ylines = False
    gl.right_labels = False
    gl.top_labels = False

    # Essas linhas configuram as linhas de grade no mapa. Elas ativam a exibição de rótulos, desativam as linhas de grade nas direções x e y, e desativam os rótulos nas bordas direita e superior do mapa.

    levels = np.linspace(0.0, 400.0, 15)
    cnplot = ax.contourf(data.longitude, data.latitude, data.sel(season=season),levels=levels, cmap='YlGnBu', extend='max')

    # Estas linhas configuram os níveis de contorno e criam um mapa de contorno usando a função contourf do Matplotlib.
    
    ax.set_title(f'{year}', fontsize=10)

    return cnplot

# Cria uma figura e uma grade de subplots
fig = plt.figure(figsize=(15, 15))
gs = gridspec.GridSpec(3, 3, hspace=0.4, wspace=0.4)

# Cria e plota os subplots para cada temporada
for i, (season, year) in enumerate(zip(primeiras_jja, range(1991, 2000))):
    row = i // 3
    col = i % 3
    create_plot(plt.subplot(gs[row, col], projection=ccrs.PlateCarree()), sclim, -44, season, year)

# Adiciona uma barra de cores
cax = plt.axes([0.2, 0.065, 0.6, 0.02])
cbar = plt.colorbar(create_plot(plt.subplot(gs[0, 0], projection=ccrs.PlateCarree()), sclim, -44, primeiras_jja[0], 1991), cax=cax, orientation='horizontal', pad=0.4)

# Salva a figura e a exibe
plt.savefig('CHIRPSMediaAnual.png', dpi=300)
plt.show()

# Abre as imagens salvas
img = Image.open("Media dos anos(CHIRPS Media Anual).png")
img1 = Image.open("CHIRPSMediaAnual.png")

# Redimensiona as imagens
img_size = img.resize((670, 1024))
img1_size = img1.resize((670, 1024))

# Cria uma nova imagem branca
img2 = Image.new("RGB", (1340, 1024), "white")

# Cola as imagens na nova imagem
img2.paste(img_size, (0, 0))
img2.paste(img1_size, (670, 0))

# Salva a nova imagem
img2.save('Medias-anual.png')
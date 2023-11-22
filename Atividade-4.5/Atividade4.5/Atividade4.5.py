import xarray as xr
import numpy as np
import geopandas as gpd
from PIL import Image

import matplotlib.pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
import matplotlib.gridspec as gridspec

estados = gpd.read_file('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/BR_UF_2022.shp')

# print(estados.head())

# # # Extraindo as dimensões
dset = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/gpcp3.2_1983-2022_sa.nc')
# print(dset)
# print(dset['time_bnds'])


# # Extrai as coordenadas de latitude e longitude do arquivo de dados
lat = dset.lat.values
lon = dset.lon.values

# Seleciona a variável 'precip' do arquivo de dados
var = dset['sat_gauge_precip']


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
levels = np.linspace(0.0, 24.0, 11)

# Inicializa um gráfico de contorno preenchido com os dados de precipitação (clim)
cnplot = ax.contourf(clim.lon, clim.lat, clim,levels=levels, cmap='YlGnBu', extend='max')

# Adiciona uma barra de cores (colorbar) na parte inferior do gráfico
cbar = plt.colorbar(cnplot, orientation='horizontal', pad=0.07, shrink=0.6)
cbar.set_label('precipitation (mm/day) \n GPCC(2002)')

# Define o título do gráfico
ax.set_title(' (matplotlib)')

# Plota os limites dos estados (supondo que sejam limites de estados geográficos)
estados.plot(ax=ax, color='none', edgecolor='black')

# Salva a figura em um arquivo chamado 's1-matplot.png' com uma resolução de 300 dpi
plt.savefig('s1-matplot.png', dpi=300)

# Exibe o gráfico na tela
plt.show()

# Define a função 'create_plot' para criar subplots no mapa
def create_plot(ax, nrow, ncol, data, tlon, season=''):
    ax.set_extent((-90.0, -30.0, -60.0, 15.0))  # Define a extensão geográfica
    ax.coastlines(resolution='110m', color='black')  # Adiciona a linha da costa
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')  # Adiciona as fronteiras

    gl = ax.gridlines(draw_labels=True)
    gl.xlines = False
    gl.ylines = False
    gl.right_labels = False
    gl.top_labels = False

    levels = np.linspace(0.0, 18.0, 11)  # Define os níveis de contorno Aqui, você está criando uma lista de 11 valores espaçados uniformemente entre 0 e 100. Esses valores representam os níveis de contorno que serão usados no gráfico.
    cnplot = ax.contourf(data.lon, data.lat, data.sel(season=season),levels = levels,  cmap='YlGnBu', extend='max')  # Cria um plot de contorno
    ax.text(tlon, -55, season, bbox=dict(facecolor='white', alpha=0.7))  # Adiciona o nome da estação

    return cnplot

# def create_plot(ax, nrow, ncol, data, tlon, season=''):

# Define uma função chamada create_plot que recebe como argumentos:
# ax: O eixo em que o gráfico será desenhado.
# nrow: O número da linha do subplot.
# ncol: O número da coluna do subplot.
# data: Os dados que serão plotados.
# tlon: A coordenada de longitude onde o texto será adicionado.
# season: O nome da estação do ano (opcional).
# ax.set_extent((-90.0, -30.0, -60.0, 15.0)):

# Define a extensão geográfica do subplot. Neste caso, limita a área geográfica a aproximadamente a região da América do Sul.
# ax.coastlines(resolution='110m', color='black'):

# Adiciona as linhas costeiras ao gráfico. resolution='110m' define a resolução da linha da costa, e color='black' define a cor como preta.
# ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black'):

# Adiciona as fronteiras (bordas) ao gráfico. linewidth=0.5 define a largura da linha e edgecolor='black' define a cor da linha como preta.
# gl = ax.gridlines(draw_labels=True):

# Adiciona as linhas de grade ao gráfico. draw_labels=True indica que rótulos serão adicionados às linhas de grade.
# gl.xlines = False

# Desativa as linhas de grade verticais.
# gl.ylines = False

# Desativa as linhas de grade horizontais.
# gl.right_labels = False

# Desativa os rótulos do lado direito.
# gl.top_labels = False

# Desativa os rótulos superiores.
# levels = np.linspace(0.0, 18.0, 11):

# Cria uma lista de 11 valores espaçados uniformemente entre 0.0 e 18.0. Esses valores representam os níveis de contorno que serão usados no gráfico.
# cnplot = ax.contourf(data.lon, data.lat, data.sel(season=season), levels=levels, cmap='YlGnBu', extend='max'):

# Cria um gráfico de contorno preenchido (contourf) utilizando os dados fornecidos (data) para a coordenada de longitude (data.lon), a coordenada de latitude (data.lat), e seleciona a estação do ano específica (data.sel(season=season)). Os níveis de contorno são definidos pela lista levels, o mapa de cores é 'YlGnBu', e extend='max' indica que a barra de cores será estendida para valores acima do valor máximo da escala de cores.
# ax.text(tlon, -55, season, bbox=dict(facecolor='white', alpha=0.7)):

# Adiciona um texto ao gráfico. O texto é a estação do ano (season) posicionado em tlon para a coordenada de longitude e -55 para a coordenada de latitude. A caixa ao redor do texto tem uma cor de fundo branca (facecolor='white') com uma opacidade de 0.7 (alpha=0.7).
# return cnplot:

# Retorna o objeto de plotagem (cnplot) que pode ser utilizado para adicionar uma barra de cores.

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
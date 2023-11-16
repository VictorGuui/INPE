import xarray as xr
import numpy as np
import geopandas as gpd
from PIL import Image

import matplotlib.pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
import matplotlib.gridspec as gridspec

estados = gpd.read_file('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/BR_UF_2022.shp')

# Carregar os dados
dset = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/arquivo_climatologia_GPCC_1991-2020.nc')

# Extraindo as dimensões de estão dentro do arquivo, que pode ser visualizada usando print(dset)
lat = dset.lat.values
lon = dset.lon.values

# Variáveis a serem plotadas que existem dentro do arquvio, que pode ser visualizada usando print(dset)
variaveis = [
    'sat_gauge_precip',
    'sat_gauge_error',
    'satellite_precip',
    'satellite_source',
    'gauge_precip',
    'probability_liquid_phase',
    'gauge_relative_weight',
]

# Configuração da projeção do mapa e sua extensão geográfica, define uma função setup_map que configura a projeção do mapa e sua extensão geográfica usando a biblioteca Cartopy.
def setup_map(ax):
    # Define a extensão geográfica do mapa
    ax.set_extent((-90.0, -30.0, -60.0, 15.0))

    # Adiciona linhas costeiras ao mapa
    ax.coastlines(resolution='110m', color='black')

    # Adiciona fronteiras políticas ao mapa
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')

    # Adiciona linhas de grade ao mapa e armazena as configurações em 'gl'
    gl = ax.gridlines(draw_labels=True)

    # Desativa as linhas de grade verticais
    gl.xlines = False

    # Desativa as linhas de grade horizontais
    gl.ylines = False

    # Desativa os rótulos do lado direito das linhas de grade
    gl.right_labels = False

    # Desativa os rótulos no topo das linhas de grade
    gl.top_labels = False


# Criação da figura
# Criação da figura com o tamanho especificado
fig = plt.figure(figsize=(15, 10))

# Cria uma grade de subplots com 2 linhas, 4 colunas e ajustes de espaçamento
gs = gridspec.GridSpec(2, 4, hspace=0.3, wspace=0.2)

# Loop sobre as variáveis
for i, var_name in enumerate(variaveis):
    # Calcula a linha e coluna do subplot com base no índice 'i'
    row = i // 4
    col = i % 4

    # Seleciona a variável específica no conjunto de dados
    var_data = dset[var_name]

    # Calcula a média anual dos dados ao longo do tempo
    clim = np.mean(var_data, axis=0)

    # Inicializa um subplot na posição 'row', 'col' com a projeção PlateCarree
    ax = plt.subplot(gs[row, col], projection=ccrs.PlateCarree())

    # Configuração do mapa usando a função 'setup_map'
    setup_map(ax)

    # Plota os dados no mapa com cores de contorno
    cnplot = ax.contourf(lon, lat, clim, cmap='YlGnBu', extend='max')

    # Adiciona uma barra de cores horizontal abaixo do subplot
    cbar = plt.colorbar(cnplot, orientation='horizontal', pad=0.07, shrink=0.6)
    cbar.set_label(f'{var_name}')

    # Mascara dos estados
    estados.plot(ax=ax, color='none', edgecolor='black')

    # Define o título do subplot com base na variável
    ax.set_title(f'{var_name} Climatology')

# Salva a figura em um arquivo chamado 'GPCCjul2005_jun2006-matplot.png' com uma resolução de 300 dpi
plt.savefig('./GPCC/GPCCjul2005_jun2006-matplot-mm.png', dpi=300)

# Exibe o gráfico na tela
plt.show()


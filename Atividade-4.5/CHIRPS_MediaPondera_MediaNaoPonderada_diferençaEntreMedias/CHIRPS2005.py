# Importação das bibliotecas necessárias
import xarray as xr    # Importa a biblioteca xarray e a apelida como "xr"
import numpy as np     # Importa a biblioteca numpy e a apelida como "np"
import geopandas as gpd # Importa a biblioteca geopandas e a apelida como "gpd"
from PIL import Image  # Importa a classe Image da biblioteca PIL (Python Imaging Library)

import matplotlib.pyplot as plt                 # Importa a biblioteca matplotlib e a apelida como "plt"
from cartopy import crs as ccrs, feature as cfeature  # Importa partes específicas da biblioteca cartopy

import matplotlib.gridspec as gridspec  # Importa a biblioteca gridspec do matplotlib

# Leitura dos estados do shapefile
estados = gpd.read_file('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/BR_UF_2022.shp')
# Lê as informações geográficas dos estados do Brasil a partir de um arquivo shapefile

# Leitura do conjunto de dados de precipitação
dset = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/arquivo_climatologia_chirps_AmericaDoSul_coord_-90,-30,-60,15_2005.nc')
# Lê um conjunto de dados climatológicos de precipitação para a América do Sul a partir de um arquivo NetCDF

# Calcular o comprimento do mês para ponderação
month_length = dset.time.dt.days_in_month
# Calcula o número de dias em cada mês

# Calcular os pesos agrupados por 'time.season'
weights = (month_length.groupby("time.season") / month_length.groupby("time.season").sum())
# Calcula os pesos para cada estação do ano com base no comprimento dos meses

# Calcular a média ponderada
dset_weighted = (dset['precip'] * weights).groupby("time.season").sum(dim="time")
# Calcula a média ponderada da precipitação para cada estação do ano

# Calcular a média não ponderada para comparação
dset_unweighted = dset['precip'].groupby("time.season").mean("time")
# Calcula a média não ponderada da precipitação para cada estação do ano

# Diferença entre média ponderada e não ponderada
dset_diff = dset_weighted - dset_unweighted
# Calcula a diferença entre a média ponderada e não ponderada da precipitação para cada estação do ano

# Configuração dos gráficos
fig, axes = plt.subplots(nrows=4, ncols=3, figsize=(14, 12))
# Configura uma matriz de subgráficos com 4 linhas, 3 colunas e tamanho total de 14x12 polegadas

levels = np.linspace(0.0, 450.0, 11)

# Loop sobre as estações para plotar os resultados
for i, season in enumerate(("DJF", "MAM", "JJA", "SON")):
    
    difLevels = np.linspace(-9.0, 9.0, 11)

    # Plotagem da média ponderada
    dset_weighted.sel(season=season).plot.pcolormesh( 
        ax=axes[i, 0], #
        cmap="YlGnBu",
        add_colorbar=True,
        extend="both",
        levels=levels,
    )
    # Plota o mapa da média ponderada da precipitação para a estação atual no primeiro subgráfico

# dset_weighted.sel(season=season): Isola os dados da variável dset_unweighted para a estação especificada pelo valor da variável season.

# .plot.pcolormesh(...): Gera um gráfico de pseudocor ou pseudocor map usando a função pcolormesh da biblioteca de plotagem associada ao conjunto de dados.

# ax=axes[i, 1]: Especifica o eixo no qual o gráfico será plotado. O valor axes[i, 1] indica a posição do subplot na matriz de eixos especificada pelo índice i na coluna 1.

# cmap="YlGnBu": Define o mapa de cores para o gráfico. Neste caso, o mapa de cores usado é "YlGnBu" (amarelo, verde, azul).

# add_colorbar=True: Adiciona uma barra de cores ao gráfico, que fornece uma referência visual para a escala de cores usada no mapa de cores.

# extend="both": Indica que a barra de cores deve ser estendida nas duas extremidades, indicando valores fora dos limites dos dados plotados.

# levels=levels: Define os níveis (intervalos) para os quais as cores do mapa de cores são mapeadas. Isso permite uma personalização mais detalhada da representação visual dos dados.

    # Plotagem da média não ponderada
    dset_unweighted.sel(season=season).plot.pcolormesh(
        ax=axes[i, 1],
        cmap="YlGnBu",
        add_colorbar=True,
        extend="both",
        levels=levels,
    )
    # Plota o mapa da média não ponderada da precipitação para a estação atual no segundo subgráfico

    # Plotagem da diferença entre ponderada e não ponderada
   # Plotagem da diferença entre ponderada e não ponderada
    dset_diff.sel(season=season).plot.pcolormesh(
    ax=axes[i, 2],
    cmap="RdBu_r",
    add_colorbar=True,
    extend="both",
)

    # Plota o mapa da diferença entre a média ponderada e não ponderada da precipitação para a estação atual no terceiro subgráfico

    axes[i, 0].set_ylabel(season)  # Adiciona o nome da estação como rótulo no eixo y no primeiro subgráfico
    axes[i, 1].set_ylabel("")  # Remove o rótulo do eixo y no segundo subgráfico
    axes[i, 2].set_ylabel("")  # Remove o rótulo do eixo y no terceiro subgráfico

    # Adiciona a máscara dos estados a todos os subgráficos
    estados.plot(ax=axes[i, 0], color='none', edgecolor='black')
    estados.plot(ax=axes[i, 1], color='none', edgecolor='black')
    estados.plot(ax=axes[i, 2], color='none', edgecolor='black')

# dset_weighted.sel(season=season).plot.pcolormesh(...): Plota o mapa da média ponderada da precipitação para a estação atual no primeiro subgráfico.

# dset_unweighted.sel(season=season).plot.pcolormesh(...): Plota o mapa da média não ponderada da precipitação para a estação atual no segundo subgráfico.

# dset_diff.sel(season=season).plot.pcolormesh(...): Plota o mapa da diferença entre a média ponderada e não ponderada da precipitação para a estação atual no terceiro subgráfico.

# o terceiro subgráfico (axes[i, 2]) mostra a diferença entre as médias ponderada e não ponderada para cada estação do ano. 

# Ajuste das configurações dos eixos
for ax in axes.flat:
    ax.axes.get_xaxis().set_ticklabels([])  # Remove os rótulos do eixo x
    ax.axes.get_yaxis().set_ticklabels([])  # Remove os rótulos do eixo y
    ax.axes.axis("tight")  # Ajusta os limites dos eixos
    ax.set_xlabel("")  # Remove o rótulo do eixo x

# Adiciona títulos aos subplots
axes[0, 0].set_title("Média Ponderada")  # Adiciona um título ao primeiro subgráfico
axes[0, 1].set_title("Média Não Ponderada")  # Adiciona um título ao segundo subgráfico
axes[0, 2].set_title("Diferença")  # Adiciona um título ao terceiro subgráfico

plt.tight_layout()  # Ajusta o layout para evitar sobreposição

plt.savefig('./CHIRPS_MediaPondera_MediaNaoPonderada_diferençaEntreMedias/CHIRPS2000_mm-mes_matplot.png', dpi=300)

# Título geral
fig.suptitle("Precipitação Média por Estação", fontsize=16, y=1.02)
# Adiciona um título geral ao gráfico

# Exibe o gráfico
plt.show()
# Mostra o gráfico na tela








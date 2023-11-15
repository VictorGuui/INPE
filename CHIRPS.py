# Importação das bibliotecas necessárias
import xarray as xr
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from cartopy import crs as ccrs, feature as cfeature
import matplotlib.gridspec as gridspec

# Leitura dos estados do shapefile
estados = gpd.read_file('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/BR_UF_2022.shp')

# Leitura do conjunto de dados de precipitação
dset = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/Chirps1AmericaDoSul1991_2020.nc')
dsetMediaAnual = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/Chirps1AmericaDoSul1991_2020MediaAnual.nc')
dsetMediaSazonal = xr.open_dataset('/mnt/c/INPE/Semana30102023-01112023/Atividade3.2.MalhaComTodasUFsDoBrasil/Chirps1AmericaDoSul1991_2020MediaSazonal.nc')

# Definindo as estações
seasons = ["DJF", "MAM", "JJA", "SON"]

season_months = {
    "DJF": [12, 1, 2],
    "MAM": [3, 4, 5],
    "JJA": [6, 7, 8],
    "SON": [9, 10, 11]
}

# Criando a grade 4x3
fig = plt.figure(figsize=(15, 12))
gs = gridspec.GridSpec(4, 3, height_ratios=[1, 1, 1, 1], width_ratios=[1, 1, 1], hspace=0.3, wspace=0.3)

# Loop para preencher as subplots
# Loop para preencher as subplots
# Loop para preencher as subplots
for i, data in enumerate([dset, dsetMediaAnual, dsetMediaSazonal]):
    for j, season in enumerate(seasons):
        print(f"Dataset variables: {list(data.variables)}")
        
        ax = plt.subplot(gs[j, i], projection=ccrs.PlateCarree())
        
        # Verificar se 'precip' está presente no conjunto de dados
        if 'precip' not in data.variables:
            print("'precip' not found in the dataset.")
            continue

        # Plotando os dados
        data_var = data['precip']
        
        # Convert season label to months
        months = season_months[season]
        
        # Select data for the specific season using the month values
        data_var_season = data_var.sel(time=data_var['time.month'].isin(months))

        # Usar pcolormesh em vez de plot para dados tridimensionais
        im = ax.pcolormesh(data_var['longitude'], data_var['latitude'], data_var_season, transform=ccrs.PlateCarree(), cmap='Blues')
        
        # Adicionar contornos dos estados
        estados.boundary.plot(ax=ax, linewidth=1, color='black')
        
        # Adicionando títulos
        ax.set_title(f"{data_var.long_name} - {season}")
        
        # Adicionando grade de latitudes e longitudes
        ax.gridlines(draw_labels=True)

# Adicionando barra de cores
cbar_ax = fig.add_axes([0.95, 0.15, 0.01, 0.7])
cbar = plt.colorbar(im, cax=cbar_ax)
cbar.set_label(data_var.units)

# Adiciona títulos aos subplots
ax[0, 0].set_title("Precip 1991/2005")
ax[0, 1].set_title("Média Anual")
ax[0, 2].set_title("Média Sazonal")

plt.show()

import utils as ut
import pandas as pd
import numpy as np
from paths import get_data_path, get_repository_path, get_visuals_path

data = pd.read_csv(get_data_path() / 'Composite_TOI_Host.csv')
ra = data['ra']
dec = data['dec']
AB = data['AB_Mag']
dist = data['sy_dist']
Teff = data['st_teff']
weight = data['Weight']
data_size = data['Data_Size']
SAT_Mag_dir = get_repository_path() / 'systematics' / 'SatMag.csv'
LIM_Mag_dir = get_repository_path() / 'systematics' / 'LimMag.csv'
CF_Map = ut.cost_function_map(ra, dec, AB, dist, Teff, weight, data_size, SAT_Mag_dir, LIM_Mag_dir, radius = 7, output_dir=None, 
                              resolution = 15, dust_map_version='bayestar2019')
# dust_map = ut.dust_map(100, dist = 500)
ut.Plot_Cartesian_Sky_Map(
    CF_Map,
    heatmap_label='Weighted ULTRASAT target yield',
    output_path=get_visuals_path() / 'ultrasat_exoplanet_field_cost.png',
)
import utils as ut
import pandas as pd
import numpy as np
data = pd.read_csv('data/Composite_TOI_Host.csv')
ra = data['ra']
dec = data['dec']
AB = data['AB_Mag']
dist = data['sy_dist']
Teff = data['st_teff']
weight = data['Weight']
data_size = data['Data_Size']
SAT_Mag_dir = '../ultrapipe/systematics/SatMag.csv'
LIM_Mag_dir = '../ultrapipe/systematics/SatMag.csv'
CF_Map = ut.cost_function_map(ra, dec, AB, dist, Teff, weight, data_size, SAT_Mag_dir, LIM_Mag_dir, radius = 7, output_dir=None, 
                              resolution = 15, dust_map_version='bayestar2019')
# dust_map = ut.dust_map(100, dist = 500)
ut.Plot_Cartesian_Sky_Map(CF_Map)
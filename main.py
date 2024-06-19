import numpy as np
import matplotlib.pyplot as plt
import sys
from landlab.components import OverlandFlow
from landlab.components import SoilInfiltrationGreenAmpt
from landlab.io import read_esri_ascii

grid_path = './Inputs/LuckyHills103_1m.asc'
outlet_node = int(14504)

grid, data = read_esri_ascii(grid_path)
grid.set_watershed_boundary_condition(node_data=data, nodata_value=-9999.0)
topo = grid.add_zeros('topographic__elevation', at='node')
soil = grid.add_zeros('soil__depth', at='node')
bedrock = grid.add_zeros('bedrock__elevation', at='node')

# Update topography
topo[:] = data
soil[:] = 1
bedrock[:] = topo-soil


# Infiltration component
hydraulic_conductivity_sediment = 5.5*10**-6
infilitration_depth = grid.add_ones("soil_water_infiltration__depth", at="node", dtype=float)
grid.add_zeros('surface_water__depth',
               at="node")
infilitration_depth *= 0.03
init_soil_moisture = 0.2
volume_fraction_coarse_fragments = 0.2

Ks = grid.add_ones('hydraulic_conductivity', at="node", dtype=float)
Ks *= hydraulic_conductivity_sediment
SI = SoilInfiltrationGreenAmpt(grid, hydraulic_conductivity=Ks,
                                     soil_type='sandy loam',
                                     initial_soil_moisture_content=init_soil_moisture ,
                                     volume_fraction_coarse_fragments=volume_fraction_coarse_fragments,
                                     )

## Overland flow component

roughness = 0.07
ld = OverlandFlow(grid, mannings_n = roughness,
                  steep_slopes=True)

# Rainfall data
rainfall_data = np.load('./Inputs/rainfall_data.npz')
rainfall_duration = rainfall_data['durations']
rainfall_rate = rainfall_data['rates']


int_index = 0
current_rainfall_duration = rainfall_duration[int_index]
ld.rainfall_intensity = rainfall_rate[int_index]
epsilon = 10**-10

elapse_dts = 0
min_dt = 30     # maximal dt for ensure stability
while elapse_dts < rainfall_duration[-1]:
    if elapse_dts >= current_rainfall_duration:  # sec

        int_index += 1
        current_rainfall_rate = rainfall_rate[int_index]
        current_rainfall_duration = rainfall_duration[int_index]
        ld.rainfall_intensity = current_rainfall_rate  # meter per sec -> For the OverlandFlow component

    ld.calc_time_step()
    dt = np.min((ld._dt,min_dt))
    ld.run_one_step(dt=dt)
    SI.run_one_step(dt=dt)
    elapse_dts += dt
    print(elapse_dts)

# Let the watershed run-out of water
ld.rainfall_intensity = epsilon
while np.max(grid.at_node['surface_water__depth']) >= 0.0001:

    ld.calc_time_step()
    dt = ld._dt
    ld.run_one_step(dt=dt)
    SI.run_one_step(dt=dt)
    elapse_dts += dt
    print(elapse_dts)



import numpy as np
from landlab.components import OverlandFlow
from landlab.components import SoilInfiltrationGreenAmpt
from landlab.io import read_esri_ascii
import matplotlib.pyplot as plt
# First, read the topo data
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


# Infiltration component init
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

## Overland flow component init
roughness = 0.07
grid.add_zeros('water_surface__slope',
               at="node")
grid.add_zeros('surface_water__depth_at_link',
               at="link")
grid.add_zeros('water_surface__elevation',
               at="node")

grid.at_link['surface_water__depth_at_link'][:] = 10 ** -8
grid.at_node['surface_water__depth'][:] = 10 ** -8

ld = OverlandFlow(grid, mannings_n = roughness,
                  steep_slopes=True)

# Read rainfall data
rainfall_data = np.load('./Inputs/rainfall_data.npz')
rainfall_duration = rainfall_data['durations']  # Seconds
rainfall_rate = rainfall_data['rates']          # Rainfall intensity [m/s]

int_index = 0   # First rainfall intensity interval
current_rainfall_duration = rainfall_duration[int_index]
ld.rainfall_intensity = rainfall_rate[int_index]
epsilon = 10**-10
threshold_of_pits = 1 # m^3

# Main loop
elapse_dts = 0  # counter of simulation time [sec]
min_dt = 30     # maximal dt [sec] to ensure stability
wh_at_node = []
while elapse_dts < rainfall_duration[-1]:
    if elapse_dts >= current_rainfall_duration:  # sec

        int_index += 1
        current_rainfall_rate = rainfall_rate[int_index]
        current_rainfall_duration = rainfall_duration[int_index]
        ld.rainfall_intensity = current_rainfall_rate  # meter per sec -> For the OverlandFlow component
    wh_at_node.append(grid.at_node['surface_water__depth'][outlet_node+1])
    ld.calc_time_step()
    dt = np.min((ld._dt,min_dt))
    ld.run_one_step(dt=dt)
    SI.run_one_step(dt=dt)
    elapse_dts += dt
    print(elapse_dts)

# Let the watershed run-out of water
ld.rainfall_intensity = epsilon # very small number
while np.sum(grid.at_node['surface_water__depth']) >= threshold_of_pits: # continue until reaching threshold of topographic 'pits' water storage

    ld.calc_time_step()
    dt = ld._dt
    ld.run_one_step(dt=dt)
    SI.run_one_step(dt=dt)
    elapse_dts += dt
    print(np.sum(grid.at_node['surface_water__depth']))
    wh_at_node.append(grid.at_node['surface_water__depth'][outlet_node+1])


plt.plot(wh_at_node),plt.show()
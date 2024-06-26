import numpy as np
from landlab.components import OverlandFlow
from landlab.components import SoilInfiltrationGreenAmpt
from landlab.io import read_esri_ascii
import matplotlib.pyplot as plt
import time

# First, read the topo data
grid_path = './Inputs/LuckyHills103_1m.asc'
outlet_node = int(14504)

grid, data = read_esri_ascii(grid_path)
grid.set_watershed_boundary_condition(node_data=data, nodata_value=-9999.0)
topo = grid.add_zeros('topographic__elevation', at='node')
soil = grid.add_zeros('soil__depth', at='node')
bedrock = grid.add_zeros('bedrock__elevation', at='node')
n_grain_sizes = 3
n_links = 4

# Create array that store the flux at link per grain size
flux_at_link_per_grain_size = np.zeros((np.size(grid.zeros(at="link")), n_grain_sizes ))

# Create array that store the flux at node for all out/in links per grain size
outlinks_fluxes_at_node = np.zeros((np.size(grid.zeros(at="node")), n_links, n_grain_sizes))
inlinks_fluxes_at_node = np.zeros((np.size(grid.zeros(at="node")), n_links, n_grain_sizes))

# Update topography
topo[:] = data
soil[:] = 1
bedrock[:] = topo-soil


# Just for the example -- apriori set constant flux at all links and all grain sizes
flux_at_link_per_grain_size[:] = 1

# Calc the topo gradient
topo_grad = grid.calc_grad_at_link('topographic__elevation')

run_times  = []
for _ in range(1000):
    t1 = time.time()
    # Find the outlinks for each node.
    outlinks_at_node = grid.link_at_node_is_downwind(topo_grad)
    inlinks_at_node = grid.link_at_node_is_upwind(topo_grad)


    # Weight flux for outlinks and inlinks at node
    fluxes_at_node = flux_at_link_per_grain_size[grid.links_at_node, :]
    outlinks_fluxes_at_node[outlinks_at_node, :] = fluxes_at_node[outlinks_at_node, :]
    inlinks_fluxes_at_node[ inlinks_at_node, :] = fluxes_at_node[inlinks_at_node, :]
    outlinks_id = grid.links_at_node[outlinks_at_node]
    tend = time.time()-t1
    run_times.append(tend)

print('Average run time: ', np.round(np.mean(run_times),5), ' Min run time: ', np.round(np.min(run_times),5), ' Max run time: ', np.round(np.max(run_times),5))





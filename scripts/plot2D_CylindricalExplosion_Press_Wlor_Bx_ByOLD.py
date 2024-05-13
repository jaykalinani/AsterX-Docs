import numpy as np
from matplotlib import pyplot as plt
import openpmd_api as io

import os
home = os.environ["HOME"]

load_str_it="it00000270"
record    = "hydrobase_rho_rl00"
component = "hydrobase_rho"

series = io.Series("/Users/jaykalinani/ET/simulations/AsterX/Cylindrical_blast/Cylindrical_blast."+load_str_it+".bp", io.Access.read_only)


print("The Series contains: ",series.iterations.size())



series.flush()
print(series)
it = series.read_iterations
str_it = str(series.iterations)

print(it)




iter_rec_comp_dict = {}
iter_rec_comp_dict[it] = {}

print(it)
for key in it.meshes:
  record = it.meshes[key]
  iter_rec_comp_dict[str_it][key] = {}
  for component in record:
    iter_rec_comp_dict[str_it][key][component] = record[component].load_chunk()

series.flush()

record    = "hydrobase_rho_rl00"
component = "hydrobase_rho"
array3D = iter_rec_comp_dict[str_it][record][component]


# Set title and labels
axplot.set_title("Rest-mass density", fontsize = 10., fontweight = "bold", color = "midnightblue")
axplot.set_xlabel("x", fontsize = 7.)
axplot.set_ylabel("y", fontsize = 7.)
axplot.tick_params(labelsize=7)
axplot.xaxis.set_major_locator(plt.MaxNLocator(5))
axplot.yaxis.set_major_locator(plt.MaxNLocator(5))



# Set up the axes for the plot and the colorbar
fig    = plt.figure(figsize = [4.75, 3.5])
axplot = fig.add_axes([0.12, 0.14, 0.75, 0.75])
axclb  = fig.add_axes([0.88, 0.14, 0.02, 0.75])

z_index = int(array3D.shape[0]/2)
x0     = np.linspace(-0.5, 0.5, array3D.shape[2])
y0     = np.linspace(-0.5, 0.5, array3D.shape[1])
image   = axplot.pcolormesh(x0, y0, array3D[z_index, :, :],
                            cmap = "magma", vmin = 0.4, vmax = 2.5)

# Set up the colorbar
axclb.tick_params(labelsize=7.0)
fig.colorbar(image, cax = axclb, extend = "neither")

# Print the current iteration
axplot.text(0.18, 0.42, "Iteration " + str_it,
  fontsize = 8., fontweight = "bold", color = "white")

fig.savefig(home+"/Desktop/CylExp_rho.pdf")

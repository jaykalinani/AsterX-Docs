import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import openpmd_api as io
import os

home = os.environ["HOME"]


dr="/Users/jaykalinani/ET/simulations/AsterX/Magnetic_loop_advection_minmod/Magnetic_loop_advection_minmod"

allruns = [
(dr, ".it00000000.bp", 000),
(dr, ".it00024576.bp", 24576),
]

ext = [-0.502, 0.502, -0.502, 0.502]

def print_data(series):
  for index in series.iterations:
    iteration = str(index)
    print("\nIteration " + iteration + ":")
    print("==============")
    i = series.iterations[index]

    for key in i.meshes:
        print("Components of record \"" + key + "\":")
        record = i.meshes[key]

        for component in record:
            print("    > " + component)  # 'component' is a string 

def load_data(series, it):
  i = series.iterations[it]
 
  # velx
  key = "hydrobasex_vel_patch00_lev00"
  component = "hydrobasex_velx"

  record = i.meshes[key]
  # Load all components of the record as data chunks
  iter_rec_comp_dict = {str(it): {key: {component: record[component].load_chunk()}}}
  # Flush the series after loading all data
  series.flush()

  # Retrieve the 3D and 2D array containing the data
  velx3D = iter_rec_comp_dict[str(it)][key][component]
  z_index = int(velx3D.shape[0] / 2)
  velx2D = velx3D[z_index, :, :]

  # vely
  key = "hydrobasex_vel_patch00_lev00"
  component = "hydrobasex_vely"

  record = i.meshes[key]
  # Load all components of the record as data chunks
  iter_rec_comp_dict = {str(it): {key: {component: record[component].load_chunk()}}}
  # Flush the series after loading all data
  series.flush()

  # Retrieve the 3D and 2D array containing the data
  vely3D = iter_rec_comp_dict[str(it)][key][component]
  z_index = int(vely3D.shape[0] / 2)
  vely2D = vely3D[z_index, :, :]

  # velz
  key = "hydrobasex_vel_patch00_lev00"
  component = "hydrobasex_velz"

  record = i.meshes[key]
  # Load all components of the record as data chunks
  iter_rec_comp_dict = {str(it): {key: {component: record[component].load_chunk()}}}
  # Flush the series after loading all data
  series.flush()

  # Retrieve the 3D and 2D array containing the data
  velz3D = iter_rec_comp_dict[str(it)][key][component]
  z_index = int(velz3D.shape[0] / 2)
  velz2D = velz3D[z_index, :, :]

  wlor2D = 1.0/np.sqrt(1.0 - velx2D*velx2D - vely2D*vely2D - velz2D* velz2D) 
  
  # Bx
  key = "hydrobasex_bvec_patch00_lev00"
  component = "hydrobasex_bvecx"

  record = i.meshes[key]
  # Load all components of the record as data chunks
  iter_rec_comp_dict = {str(it): {key: {component: record[component].load_chunk()}}}
  # Flush the series after loading all data
  series.flush()

  # Retrieve the 3D and 2D array containing the data
  Bx3D = iter_rec_comp_dict[str(it)][key][component]
  z_index = int(Bx3D.shape[0] / 2)
  Bx2D = Bx3D[z_index, :, :]

  #By
  key = "hydrobasex_bvec_patch00_lev00"
  component = "hydrobasex_bvecy"

  record = i.meshes[key]
  # Load all components of the record as data chunks
  iter_rec_comp_dict = {str(it): {key: {component: record[component].load_chunk()}}}
  # Flush the series after loading all data
  series.flush()

  # Retrieve the 3D and 2D array containing the data
  By3D = iter_rec_comp_dict[str(it)][key][component]
  z_index = int(By3D.shape[0] / 2)
  By2D = By3D[z_index, :, :]

  # Bz
  key = "hydrobasex_bvec_patch00_lev00"
  component = "hydrobasex_bvecz"

  record = i.meshes[key]
  # Load all components of the record as data chunks
  iter_rec_comp_dict = {str(it): {key: {component: record[component].load_chunk()}}}
  # Flush the series after loading all data
  series.flush()

  # Retrieve the 3D and 2D array containing the data
  Bz3D = iter_rec_comp_dict[str(it)][key][component]
  z_index = int(Bz3D.shape[0] / 2)
  Bz2D = Bz3D[z_index, :, :]

  # b2small
  bst2D = wlor2D * (Bx2D*velx2D + By2D*vely2D + Bz2D*velz2D)
  B2big2D = Bx2D*Bx2D + By2D*By2D + Bz2D*Bz2D

  b2small2D = (B2big2D + bst2D*bst2D) / (wlor2D*wlor2D)
 
  return Bx2D, b2small2D
  
def prepare_grid():
  fig= plt.figure(tight_layout=True)
#  fig   = plt.figure(figsize=[4.5,3.5])
#  grid  = ImageGrid(fig, 111,
#          nrows_ncols = (1,1),
#          direction = "row",
#          axes_pad=(1.3, 0.3),
#          share_all=False,
#          aspect=True,
#          label_mode = "L",
#          cbar_size="3%",
#          cbar_pad="1.5%",
#          cbar_location = "right",
#          cbar_mode="each")
 
  grid  = ImageGrid(fig, 111, nrows_ncols = (2,2), share_all=False, aspect=True, axes_pad=0.12)
  return grid

def plot_Bx(data, ax):
  im = ax.imshow(data, vmin=-1e-3, vmax=1e-3, interpolation='bilinear',
                   cmap='Greys', extent=ext, origin='lower')

  ax.set_ylabel(r'$y$', fontsize=6.5)

  minorLocator = ticker.MultipleLocator(0.1)
  majorLocator = ticker.MultipleLocator(0.4)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

  ax.set_xlabel(r'$x$', fontsize=6.5)

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)

  x = [-0.08,0.498,0.498,-0.08]
  y= [0.38,0.38,0.498,0.498]
  ax.fill(x, y, "w")

  cbaxes = ax.inset_axes(
        [-0.02, 0.437, 0.465, 0.051], transform=ax.transData)
  cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal') #, ticks=[-1, 0, 1])
  cb.ax.set_xticklabels([r'$-10^{-3}$', r'$0$', r'$10^{-3}$'])
  cb.ax.tick_params(axis='x', direction='in', length=7, width=0.4,
    pad=2.5)
  for item in (cb.ax.get_xticklabels()):
      item.set_fontsize(5.25)

  timestr = r'$\mathrm{B^x}$'
  #timestr = r'$\mathrm{Press}$'
  cb.ax.text(
   -0.003 , 0.45, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5,
    zorder=2)

def plot_b2small(data, ax):
  im = ax.imshow(data/2.0, vmin=1e-7, vmax=1e-6, interpolation='bilinear',
                   cmap='Greys', extent=ext, origin='lower')
  #ccc
  #ax.contour(data/2, np.arange(1e-7, 1e-6, 0.5e-7), colors='k', axes=ax, linewidths=0.4, linestyles = 'solid')
  ax.contour(data/2, np.arange(0.1*np.max(data/2.0), np.max(data/2.0),0.1*np.max(data/2.0) ), 
             colors='k', axes=ax, linewidths=0.4, linestyles = 'solid',
             extent=ext)
  ax.set_ylabel(r'$y$', fontsize=6.5)

  minorLocator = ticker.MultipleLocator(0.1)
  majorLocator = ticker.MultipleLocator(0.4)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

  ax.set_xlabel(r'$x$', fontsize=6.5)

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)

  x = [0,0.5,0.5,0]
  y= [0.433,0.433,0.5,0.5]
  x = [-0.08,0.5,0.5,-0.08]
  y= [0.38,0.38,0.5,0.5]
#  ax.fill(x, y, "w")
  cbaxes = ax.inset_axes(
        [0.0, 0.437, 0.465, 0.051], transform=ax.transData)
  cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks=[0.1e-6, 0.33e-6, 0.667e-6, 1e-6])

  cb.ax.set_xticklabels([r'$0.1$', r'$0.33$', r'$0.667$',  r'$1.0$'])
  cb.ax.tick_params(axis='x', direction='in', length=7.5, width=0.4,
    pad=2.5)
  for item in (cb.ax.get_xticklabels()):
      item.set_fontsize(5.25)

  timestr = r'$\mathrm{Magnetic \ pressure}$'
  cb.ax.text(
   -8.5e-7 , 0.45, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5,
    zorder=2)

  timestr2 = r'$\times 10^{-6}$'
  cb.ax.text(
   -8.5e-7 , -0.45, timestr2, color='k', horizontalalignment='left', verticalalignment='top',size=6.5,
    zorder=2)
def main(runs):
  grid = prepare_grid()
  i=0
  for (dr, str_it, it) in runs:
    series = io.Series(dr+str_it, io.Access.read_only)
    print_data(series)
    print("\nBegin loading data..")
    Bx, b2small = load_data(series, it)
    plot_Bx(Bx, grid[i])
 #   plot_Bx(Bx, grid[i+1])
    plot_b2small(b2small, grid[i+1])
    i+=2
  print("\nDone plotting!")
  plt.savefig(home+"/Desktop/2D_MagLoopAdv.pdf", bbox_inches = 'tight', pad_inches = 0.03)
    
main(allruns)

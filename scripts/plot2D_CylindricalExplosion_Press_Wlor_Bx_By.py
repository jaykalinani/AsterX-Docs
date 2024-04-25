import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import openpmd_api as io
import os

home = os.environ["HOME"]


dr="/Users/jaykalinani/ET/simulations/AsterX/Cylindrical_blast/Cylindrical_blast"
allruns = [(dr, ".it00000270.bp", 270)]

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
  
  # (1) Pressure
  key = "hydrobasex_press_lev00"
  component = "hydrobasex_press"

  record = i.meshes[key]
  # Load all components of the record as data chunks
  iter_rec_comp_dict = {str(it): {key: {component: record[component].load_chunk()}}}
  # Flush the series after loading all data
  series.flush()

  # Retrieve the 3D and 2D array containing the data
  press3D = iter_rec_comp_dict[str(it)][key][component]
  z_index = int(press3D.shape[0] / 2)
  press2D = press3D[z_index, :, :]

  # (2) Wlor
  # velx
  key = "hydrobasex_vel_lev00"
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
  key = "hydrobasex_vel_lev00"
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
  key = "hydrobasex_vel_lev00"
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

  
  # (3) Bx
  key = "hydrobasex_bvec_lev00"
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

  # (4) By
  key = "hydrobasex_bvec_lev00"
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

  return press2D, wlor2D, Bx2D, By2D
  
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

#def plot_var(ax, var, vmin, vmax, cmap):
#  ext = [-6.01, 6.01, -6.01, 6.01]
#  im = ax.imshow(np.log10(var), vmin=vmin, vmax=vmax, interpolation='bilinear',
#                   cmap=cmap, extent=ext, origin='lower')
#  return im

def plot_press(data, ax):
  
  ext = [-6.01, 6.01, -6.01, 6.01]
  im = ax.imshow(np.log10(data), vmin=-5, vmax=-1, interpolation='bilinear',
                   cmap='Greys', extent=ext, origin='lower')
  ax.set_ylabel(r'$y$', fontsize=6.5)

  minorLocator = ticker.MultipleLocator(0.6)
  majorLocator = ticker.MultipleLocator(2.4)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)  

  x = [-0.71,5.95,5.95,-0.71]
  y= [4.40,4.40,5.9,5.9]

  ax.fill(x, y, "w")

  cbaxes = ax.inset_axes(
        [0.1, 5.1, 5.4, 0.6], transform=ax.transData)
  cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks=[-5, -3, -1])
  cb.ax.set_xticklabels([r'$10^{-5}$', r'$10^{-3}$', r'$10^{-1}$'])
  cb.ax.tick_params(axis='x', direction='in', length=7, width=0.4, 
    pad=2.5)
  for item in (cb.ax.get_xticklabels()):
      item.set_fontsize(5.25)

  timestr = r'$\mathrm{Pressure}$'
  cbaxes.text(-9.4 , 0.55, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5)

def plot_wlor_Blines(wlor, Bx, By, ax):
  ext = [-6.01, 6.01, -6.01, 6.01]
  im = ax.imshow(wlor, vmin=1, vmax=4.57, interpolation='bilinear',
                   cmap='Greys', extent=ext, origin='lower')
 
  x = np.linspace(-6,6,201)
  y = np.linspace(-6,6,201)
  xv, yv = np.meshgrid(x, y)
  
  #Az = Bx*yv - By*xv
  #ax.contour(Az, np.arange(-1.0, 1.0, 0.04), colors='k', axes=ax, linewidths=0.4, linestyles = 'solid') 
  
  #For streamlines
  #color = 2 * np.log(np.hypot(Bx, By))
  ax.streamplot(x, y, Bx, By, linewidth=0.4, color = 'k', # color=color, cmap=plt.cm.Greys,
              density=0.5, arrowstyle='-', arrowsize=0.3, 
              zorder=1,
              broken_streamlines=False)
              #alpha=0.8) # broken_streamlines=False)
  
  #For quiver
#  ax.quiver(x, y, Bx, By, 
#         units='xy', 
#         scale=0.8, zorder=1, color='k',
#         width=0.01, 
#         headwidth=3., 
#         headlength=4.)

  ax.set_ylabel(r'$y$', fontsize=6.5)

  minorLocator = ticker.MultipleLocator(0.6)
  majorLocator = ticker.MultipleLocator(2.4)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)

  x = [-0.71,5.95,5.95,-0.71]
  y= [4.40,4.40,5.9,5.9]

  x1 = [-5.9,-1.9,-1.9,-5.9]
  y1= [4.8,4.8,5.4,5.4]

  x2 = [-5.9,-0.72,-0.72,-5.9]
  y2= [4.20,4.20,4.82,4.82]
  
  ax.fill(x, y, "w")
  ax.fill(x1, y1, "w")
  ax.fill(x2, y2, "w")

  cbaxes = ax.inset_axes(
        [0.1, 5.1, 5.4, 0.6], transform=ax.transData)
  cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks=[1, 2.19, 3.38, 4.57])
  cb.ax.set_xticklabels([r'$1.00$', r'$2.19$', r'$3.38$',  r'$4.57$'])
  cb.ax.tick_params(axis='x', direction='in', length=7, width=0.4,
    pad=2.5)
  for item in (cb.ax.get_xticklabels()):
      item.set_fontsize(5.25)

  timestr = r'$\mathrm{Lorentz \ factor}$'
  cbaxes.text(-2.9 , 0.55, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5)
  timestr = r'$\mathrm{Magnetic \ fieldlines}$'
  cbaxes.text(-2.9 , -0.55, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5)

def plot_Bx(data, ax):
  ext = [-6.01, 6.01, -6.01, 6.01]
  im = ax.imshow(data, vmin=0, vmax=0.36, interpolation='bilinear',
                   cmap='Greys', extent=ext, origin='lower')
  
  ax.set_ylabel(r'$y$', fontsize=6.5)

  minorLocator = ticker.MultipleLocator(0.6)
  majorLocator = ticker.MultipleLocator(2.4)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

  ax.set_xlabel(r'$x$', fontsize=6.5)

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)

  x = [-0.71,5.95,5.95,-0.71]
  y= [4.40,4.40,5.9,5.9]
  ax.fill(x, y, "w")

  cbaxes = ax.inset_axes(
        [0.1, 5.1, 5.4, 0.6], transform=ax.transData)
  cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks=[0, 0.12, 0.24,0.36])
  cb.ax.set_xticklabels([r'$0.00$', r'$0.12$', r'$0.24$',  r'$0.36$'])
  cb.ax.tick_params(axis='x', direction='in', length=7, width=0.4,
    pad=2.5)
  for item in (cb.ax.get_xticklabels()):
      item.set_fontsize(5.25)

  timestr = r'$\mathrm{B^x}$'
  #timestr = r'$\mathrm{Press}$'
  cb.ax.text(
   -0.395 , 0.55, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5,
    zorder=2)


def plot_By(data, ax):
  ext = [-6.01, 6.01, -6.01, 6.01]
  im = ax.imshow(data, vmin=-0.18, vmax=0.18, interpolation='bilinear',
                   cmap='Greys', extent=ext, origin='lower')
 
  ax.set_xlabel(r'$x$', fontsize=6.5)
  ax.set_ylabel(r'$y$', fontsize=6.5)

  minorLocator = ticker.MultipleLocator(0.6)
  majorLocator = ticker.MultipleLocator(2.4)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)

  x = [-0.71,5.95,5.95,-0.71]
  y= [4.40,4.40,5.9,5.9]
  ax.fill(x, y, "w")

  cbaxes = ax.inset_axes(
        [0.1, 5.1, 5.4, 0.6], transform=ax.transData)
  cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks=[-0.18, -0.06, 0.06, 0.18])
  cb.ax.set_xticklabels([r'$-0.18$', r'$-0.06$', r'$0.06$',  r'$0.18$'])
  cb.ax.tick_params(axis='x', direction='in', length=7, width=0.4,
    pad=2.5)
  for item in (cb.ax.get_xticklabels()):
      item.set_fontsize(5.25)

  timestr = r'$\mathrm{B^y}$'
  cbaxes.text(-0.575 , 0.55, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5)



def main(runs):
  for (dr, str_it, it) in runs:
    series = io.Series(dr+str_it, io.Access.read_only)
    print_data(series)
    grid = prepare_grid()
    print("\nBegin loading data..")
    press, wlor, Bx, By = load_data(series, it)
    plot_press(press, grid[0])
    plot_wlor_Blines(wlor, Bx, By, grid[1])
    plot_Bx(Bx, grid[2])
    plot_By(By, grid[3])
    print("\nDone plotting!")
    plt.savefig(home+"/Desktop/CylBW2D.pdf", bbox_inches = 'tight', pad_inches = 0.03)

main(allruns)

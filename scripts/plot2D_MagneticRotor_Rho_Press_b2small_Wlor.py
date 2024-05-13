import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import openpmd_api as io
import os

home = os.environ["HOME"]


dr="/Users/jaykalinani/ET/simulations/AsterX/Magnetic_rotor/Magnetic_rotor"

allruns = [(dr, ".it00000640.bp", 640)]

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
  
  # (1) Rho
  key = "hydrobasex_rho_lev00"
  component = "hydrobasex_rho"

  record = i.meshes[key]
  # Load all components of the record as data chunks
  iter_rec_comp_dict = {str(it): {key: {component: record[component].load_chunk()}}}
  # Flush the series after loading all data
  series.flush()

  # Retrieve the 3D, 2D and 1D arrays containing the data
  rho3D = iter_rec_comp_dict[str(it)][key][component]
  z_index = int(rho3D.shape[0] / 2) # on xy plane

  rho2D = rho3D[z_index, :, :]

  # (2) Pressure
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

  # (3) Wlor
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

  # (4) b2/2
  # Bx
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

  #By
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

  # Bz
  key = "hydrobasex_bvec_lev00"
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
 
  return rho2D, press2D, b2small2D, wlor2D, Bx2D, By2D
  
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


def plot_rho(data, ax):

  ext = [0.01, 0.99, 0.01, 0.99]
  im = ax.imshow((data), vmin=0.2, vmax=8.2, interpolation='bilinear',
                   cmap='Greys', extent=ext, origin='lower')
  ax.set_ylabel(r'$y$', fontsize=6.5)

  minorLocator = ticker.MultipleLocator(0.04)
  majorLocator = ticker.MultipleLocator(0.2)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)

  x = [-0.71/6.,5.95/6.,5.95/6.,-0.71/6.]
  y= [4.40/6.,4.40/6.,5.9/6.,5.9/6.]

  x = [0.455,0.99,0.99,0.455]
  y= [0.867,0.867,0.99,0.99]
  ax.fill(x, y, "w")

  cbaxes = ax.inset_axes(
        [0.5, 0.92, 0.45, 0.05], transform=ax.transData)
  cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks=[0.2, 2.73, 5.46, 8.20])
  cb.ax.set_xticklabels([r'$0.20$', r'$2.73$', r'$5.46$', r'$8.20$'])
  cb.ax.tick_params(axis='x', direction='in', length=7.5, width=0.2,
    pad=2.5)
  for item in (cb.ax.get_xticklabels()):
      item.set_fontsize(5.25)

  timestr = r'$\mathrm{Density}$'
  cbaxes.text(-8.25 , 0.55, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5)



def plot_press(data, ax):
  
  ext = [0.01, 0.99, 0.01, 0.99]
  im = ax.imshow((data), vmin=0.05, vmax=3.88, interpolation='bilinear',
                   cmap='Greys', extent=ext, origin='lower')
  ax.set_ylabel(r'$y$', fontsize=6.5)

  minorLocator = ticker.MultipleLocator(0.04)
  majorLocator = ticker.MultipleLocator(0.2)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)  

  x = [-0.71/6.,5.95/6.,5.95/6.,-0.71/6.]
  y= [4.40/6.,4.40/6.,5.9/6.,5.9/6.]
  x = [0.455,0.99,0.99,0.455]
  y= [0.867,0.867,0.99,0.99]
  ax.fill(x, y, "w")

  cbaxes = ax.inset_axes(
        [0.5, 0.92, 0.45, 0.05], transform=ax.transData)
  cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks=[0.05, 1.28, 2.56, 3.88])
  cb.ax.set_xticklabels([r'$0.05$', r'$1.28$', r'$2.56$', r'$3.88$'])
  cb.ax.tick_params(axis='x', direction='in', length=7.5, width=0.4, 
    pad=2.5)
  for item in (cb.ax.get_xticklabels()):
      item.set_fontsize(5.25)

  timestr = r'$\mathrm{Pressure}$'
  cbaxes.text(-4 , 0.55, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5)

def plot_b2small(data, ax):
  ext = [0.01, 0.99, 0.01, 0.99]
  im = ax.imshow(data/2.0, vmin=0, vmax=2.43, interpolation='bilinear',
                   cmap='Greys', extent=ext, origin='lower')

  ax.set_ylabel(r'$y$', fontsize=6.5)

  minorLocator = ticker.MultipleLocator(0.04)
  majorLocator = ticker.MultipleLocator(0.2)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

  ax.set_xlabel(r'$x$', fontsize=6.5)

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)

  x = [0.455,0.99,0.99,0.455]
  y= [0.867,0.867,0.99,0.99]
  ax.fill(x, y, "w")

  cbaxes = ax.inset_axes(
        [0.5, 0.92, 0.45, 0.05], transform=ax.transData)
  cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks=[0.00, 0.81, 1.62, 2.43])
  cb.ax.set_xticklabels([r'$0.00$', r'$0.81$', r'$1.62$',  r'$2.43$'])
  cb.ax.tick_params(axis='x', direction='in', length=7.5, width=0.4,
    pad=2.5)
  for item in (cb.ax.get_xticklabels()):
      item.set_fontsize(5.25)

  timestr = r'$\mathrm{Magnetic \ pressure}$'
  #timestr = r'$\mathrm{Press}$'
  cb.ax.text(
   -2.6 , 0.55, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5,
    zorder=2)

def plot_wlor_Blines(wlor, Bx, By, ax):
  ext = [0.01, 0.99, 0.01, 0.99]
  im = ax.imshow(wlor, vmin=1, vmax=1.79, interpolation='bilinear',
                   cmap='Greys', extent=ext, origin='lower')
 
  x = np.linspace(0.01,0.99,401)
  y = np.linspace(0.01,0.99,401)
  xv, yv = np.meshgrid(x, y)
  
  #Az = Bx*yv - By*xv
  #ax.contour(Az, np.arange(-1.0, 1.0, 0.04), colors='k', axes=ax, linewidths=0.4, linestyles = 'solid') 
  
  #For streamlines
  #color = 2 * np.log(np.hypot(Bx, By))
  ax.streamplot(x, y, Bx, By, linewidth=0.4, color = 'k', # color=color, cmap=plt.cm.Greys,
              density=0.6, arrowstyle='-', arrowsize=0.3, 
              zorder=1, 
              broken_streamlines=False,
              )
              #alpha=0.8)broken_streamlines=False)
  
  #For quiver
#  ax.quiver(x, y, Bx, By, 
#         units='xy', 
#         scale=0.8, zorder=1, color='k',
#         width=0.01, 
#         headwidth=3., 
#         headlength=4.)

  ax.set_ylabel(r'$y$', fontsize=6.5)
  ax.set_xlabel(r'$x$', fontsize=6.5)
  minorLocator = ticker.MultipleLocator(0.04)
  majorLocator = ticker.MultipleLocator(0.2)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)

  x = [-0.71,5.95,5.95,-0.71]
  y= [4.40,4.40,5.9,5.9]

  x = [0.455,0.99,0.99,0.455]
  y= [0.867,0.867,0.99,0.99]
  #trials 
  x1 = [0.01,0.33,0.33,0.01]
  y1= [0.92,0.92,0.94,0.94]

  x2 = [-5.9,-0.72,-0.72,-5.9]
  y2= [4.20,4.20,4.82,4.82]
  
  x2 = [0.01,0.43,0.43,0.01]
  y2= [0.845,0.845,0.895,0.895]

  ax.fill(x, y, "w")
  ax.fill(x1, y1, "w")
  ax.fill(x2, y2, "w")

  cbaxes = ax.inset_axes(
        [0.5, 0.92, 0.45, 0.05], transform=ax.transData)
  cb = plt.colorbar(im, cax=cbaxes, orientation='horizontal', ticks=[1, 1.26, 1.52, 1.79])
  cb.ax.set_xticklabels([r'$1.00$', r'$1.26$', r'$1.52$',  r'$1.79$'])
  cb.ax.tick_params(axis='x', direction='in', length=7.5, width=0.4,
    pad=2.5)
  for item in (cb.ax.get_xticklabels()):
      item.set_fontsize(5.25)

  timestr = r'$\mathrm{Lorentz \ factor}$'
  cbaxes.text(0.15 , 0.55, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5)
  timestr = r'$\mathrm{Magnetic \ fieldlines}$'
  cbaxes.text(0.15 , -0.55, timestr, color='k', horizontalalignment='left', verticalalignment='top',size=6.5)

def main(runs):
  for (dr, str_it, it) in runs:
    series = io.Series(dr+str_it, io.Access.read_only)
    print_data(series)
    grid = prepare_grid()
    print("\nBegin loading data..")
    rho, press, b2small, wlor, Bx, By = load_data(series, it)
    plot_rho(rho, grid[0])
    plot_press(press, grid[1])
    plot_b2small(b2small, grid[2])
    plot_wlor_Blines(wlor, Bx, By, grid[3])
    #plot_Bx(Bx, grid[2])
    #plot_By(By, grid[3])
    print("\nDone plotting!")
    plt.savefig(home+"/Desktop/MagRot2D.pdf", bbox_inches = 'tight', pad_inches = 0.03)

main(allruns)

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
  y_index = int(rho2D.shape[0] / 2)
  x_index = int(rho2D.shape[1] / 2)
  rho_x = rho2D[y_index, :-1]
  rho_y = rho2D[:-1, x_index]
  
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

  y_index = int(press2D.shape[0] / 2)
  x_index = int(press2D.shape[1] / 2)
  press_x = press2D[y_index, :-1]
  press_y = press2D[:-1, x_index]


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
  y_index = int(wlor2D.shape[0] / 2)
  x_index = int(wlor2D.shape[1] / 2)
  wlor_x = wlor2D[y_index, :-1]
  wlor_y = wlor2D[:-1, x_index]
  
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
  y_index = int(b2small2D.shape[0] / 2)
  x_index = int(b2small2D.shape[1] / 2)
  b2small_x = b2small2D[y_index, :-1]
  b2small_y = b2small2D[:-1, x_index]

  return rho_x, rho_y, wlor_x, wlor_y, press_x, press_y, b2small_x, b2small_y
  
def prepare_grid():
 # fig= plt.figure() #tight_layout=True, sharex=True, sharey=True)
  fig, axs = plt.subplots(4, 2, 
            figsize=(5,4.9), 
            sharex=True,
            sharey = 'row', 
            #tight_layout=True,
            layout="constrained",
            subplot_kw=dict(box_aspect=0.47),)

#  grid  = ImageGrid(fig, 111, 
#          nrows_ncols = (4,2), 
#          share_all=False, aspect=True, axes_pad=0) #, axes_pad=0.12)
  return axs

def plot_rho(rho, ax, lb):
  
  x = np.linspace(-5.9,5.9,200)
  
  ax.plot(x, rho*1.0e4, '-', color= 'k', lw=0.8, alpha=0.9)
  if lb: 
     ax.set_ylabel(r'$\rho \times 10^{4}$', fontsize=6.5)
#  else: 
#     ax.set_yticks([])
  ax.set_ylim(0, 8)
  ax.set_xlim(-5.9, 5.9)
  minorLocator = ticker.MultipleLocator(0.5)
  majorLocator = ticker.MultipleLocator(3.0)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  minorLocator = ticker.MultipleLocator(0.5)
  majorLocator = ticker.MultipleLocator(2)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

 # ax.set_aspect(1)
  ax.grid(alpha=0.4)
  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)  

def plot_wlor(wlor, ax, lb):

  #ax.set_aspect(1)
  x = np.linspace(-5.9,5.9,200)

  ax.plot(x, wlor, '-', color= 'k', lw=0.8, alpha=0.9)
  if lb: 
     ax.set_ylabel(r'$W$', fontsize=6.5)
  ax.set_ylim(0.5, 5)
  ax.set_xlim(-5.9, 5.9)
  minorLocator = ticker.MultipleLocator(0.5)
  majorLocator = ticker.MultipleLocator(2.0)
  ax.xaxis.set_minor_locator(minorLocator)
  ax.xaxis.set_major_locator(majorLocator)
  minorLocator = ticker.MultipleLocator(0.5)
  majorLocator = ticker.MultipleLocator(1)
  ax.yaxis.set_minor_locator(minorLocator)
  ax.yaxis.set_major_locator(majorLocator)

  ax.grid(alpha=0.4)
  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)


def plot_press(press, ax, lb):
  
  x = np.linspace(-5.9,5.9,200)
  
  ax.semilogy(x, press, '-', color= 'k', lw=0.8, alpha=0.9)
  if lb:
     ax.set_ylabel(r'$ P $', fontsize=6.5)
  ax.set_ylim(9.9e-6,8e-2)
  ax.set_xlim(-5.9, 5.9)
  ax.grid(alpha=0.4)
 # minorLocator = ticker.MultipleLocator(0.5)
 # majorLocator = ticker.MultipleLocator(3.0)
 # ax.xaxis.set_minor_locator(minorLocator)
 # ax.xaxis.set_major_locator(majorLocator)
  minorLocator = ticker.MultipleLocator(5e-3)
  majorLocator = ticker.MultipleLocator(1e-2)
  #ax.yaxis.set_minor_locator(minorLocator)
   
  ax.yaxis.set_major_locator(ticker.AsinhLocator(1e-1, 101))
  #ax.yaxis.set_major_locator(ticker.LogLocator(subs='auto'))

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)

def plot_b2small(b2small, ax, lb):
  
  x = np.linspace(-5.9,5.9,200)
  
  ax.semilogy(x, b2small/2.0, '-', color= 'k', lw=0.8, alpha=0.9)
  if lb:
     ax.set_ylabel(r'$ b^2/2 $', fontsize=6.5)
     ax.set_xlabel(r'$x$', fontsize=6.5)
  else:
     ax.set_xlabel(r'$y$', fontsize=6.5)
  ax.set_ylim(9.9e-6,8e-2)
  ax.set_xlim(-5.9, 5.9)
  ax.grid(alpha=0.4)
  ax.yaxis.set_major_locator(ticker.AsinhLocator(1e-1, 101))
 # minorLocator = ticker.MultipleLocator(0.5)
 # majorLocator = ticker.MultipleLocator(3.0)
 # ax.xaxis.set_minor_locator(minorLocator)
 # ax.xaxis.set_major_locator(majorLocator)
 # minorLocator = ticker.MultipleLocator(0.5)
 # majorLocator = ticker.MultipleLocator(2)
 # ax.yaxis.set_minor_locator(minorLocator)
 # ax.yaxis.set_major_locator(majorLocator)

  for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(6.5)

def main(runs):
  for (dr, str_it, it) in runs:
    series = io.Series(dr+str_it, io.Access.read_only)
    print_data(series)
    axs = prepare_grid()
    print("\n Begin loading data..")
    rho_x, rho_y, wlor_x, wlor_y, press_x, press_y, b2small_x, b2small_y = load_data(series, it)
    
    plot_rho(rho_x, axs[0,0], 1)
    plot_rho(rho_y, axs[0,1], 0)
    plot_wlor(wlor_x, axs[1,0], 1)
    plot_wlor(wlor_y, axs[1,1], 0)
    plot_press(press_x, axs[2,0], 1)
    plot_press(press_y, axs[2,1], 0)
    plot_b2small(b2small_x, axs[3,0], 1)
    plot_b2small(b2small_y, axs[3,1], 0)

    print("\nDone plotting!")
    plt.savefig(home+"/Desktop/1D_CylBW.pdf", bbox_inches = 'tight', pad_inches = 0.03)

main(allruns)

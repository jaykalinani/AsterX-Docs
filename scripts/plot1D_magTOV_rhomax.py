import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import AxesGrid, ImageGrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import os
home = os.environ["HOME"]

dr_dx05 = "/Users/jaykalinani/Doc/PhD_Work/13_CarpetX/TOV_data_AsterX/magTOV_Z4c_AMR_dx05/"
dr_dx025 = "/Users/jaykalinani/Doc/PhD_Work/13_CarpetX/TOV_data_AsterX/magTOV_Z4c_AMR_dx025/"
dr_dx0125 = "/Users/jaykalinani/Doc/PhD_Work/13_CarpetX/TOV_data_AsterX/magTOV_Z4c_AMR_dx0125/"


rho_dx05 = np.genfromtxt(dr_dx05+"paper-hydrobasex-rho.tsv")
rho_dx05_t = rho_dx05[:,1]
rho_dx05_max = rho_dx05[:,3]
print(rho_dx05_max.shape)
lbl_dx05 = r'$\Delta x_f=0.5$'

rho_dx025 = np.genfromtxt(dr_dx025+"paper-hydrobasex-rho.tsv")
rho_dx025_t = rho_dx025[:,1]
rho_dx025_max = rho_dx025[:,3]
print(rho_dx025_max.shape)
lbl_dx025 = r'$\Delta x_f=0.25$'

rho_dx0125 = np.genfromtxt(dr_dx0125+"paper-hydrobasex-rho.tsv")
rho_dx0125_t = rho_dx0125[:,1]
rho_dx0125_max = rho_dx0125[:,3]
print(rho_dx0125_max.shape)
lbl_dx0125 = r'$\Delta x_f=0.125$'


#For convergence check
f_dx05 = interp1d(rho_dx05_t, rho_dx05_max, kind='linear')
f_dx025 = interp1d(rho_dx025_t, rho_dx025_max, kind='linear')
f_dx0125 = interp1d(rho_dx0125_t, rho_dx0125_max, kind='linear')

rho_dx05_t_new = np.linspace(0, 1670, 840)    # New x values for interpolation
rho_dx05_max_new = f_dx05(rho_dx05_t_new)

rho_dx025_t_new = np.linspace(0, 1670, 840)    # New x values for interpolation
rho_dx025_max_new = f_dx025(rho_dx025_t_new)

rho_dx0125_t_new = np.linspace(0, 1670, 840)    # New x values for interpolation
rho_dx0125_max_new = f_dx0125(rho_dx0125_t_new)

conv_lbl_dx05 = r'$|\Delta x_{f,0.25} - \Delta x_{f,0.5}|$'
conv_lbl_dx0125 = r'$|\Delta x_{f,0.125} - \Delta x_{f,0.25}|\times 2^{2.7}$'

lwd=1.3
mkr2="x"
mkrs = 3
mke=5
fs=14.5

#fig = plt.figure(figsize=(7.5,6),tight_layout=True)
fig, axs = plt.subplots(2, 1,
            figsize=(7,7.2),
            sharex='col',
            sharey = 'row',
            #tight_layout=True,
            layout="constrained",
            height_ratios = [1.8,1],
          #  width_ratios = [1,1],
          #  subplot_kw=dict(box_aspect=0.47),
            )

#plt.plot(rhos_t/203.02, (rhos_max/rhos_max[0]), lw=1.5, linestyle='-', markersize=2, marker="o", color='k', label=lbl)

ms_cu = 2.0301291587112965e+02


ax = axs[0]

#plt.plot(rho_dx05_t/ms_cu, (rho_dx05_max/rho_dx05_max[0]), lw=1.5, linestyle='-', markersize=2, marker="x", color='darkblue', label=lbl_dx05, alpha=0.7)
#plt.plot(rho_dx025_t/ms_cu, (rho_dx025_max/rho_dx025_max[0]), lw=1.5, linestyle='-', markersize=2, marker="x", color='darkgreen', label=lbl_dx025, alpha=0.7)
#plt.plot(rho_dx0125_t/ms_cu, (rho_dx0125_max/rho_dx0125_max[0]), lw=1.5, linestyle='-', markersize=2, marker="x", color='darkred', label=lbl_dx0125, alpha=0.7)

ax.plot(rho_dx05_t_new/ms_cu, (rho_dx05_max_new/rho_dx05_max_new[0]), lw=lwd, linestyle='-', markersize=2, marker="x", markevery=mke, color='tab:blue', label=lbl_dx05)#, alpha=0.6)
ax.plot(rho_dx025_t_new/ms_cu, (rho_dx025_max_new/rho_dx025_max_new[0]), lw=lwd, linestyle='-', markersize=2, marker="x", markevery=mke, color='tab:green', label=lbl_dx025) #, alpha=0.6)
ax.plot(rho_dx0125_t_new/ms_cu, (rho_dx0125_max_new/rho_dx0125_max_new[0]), lw=lwd, linestyle='-', markersize=2, marker="x", markevery=mke, color='tab:red', label=lbl_dx0125) #, alpha=0.6)

#ax.set_xlabel('t [ms]', fontsize = fs)
ax.set_ylabel(r'$\rho_{\rm c} / \rho_{\rm c, 0}$', fontsize = fs)
#ax = plt.axes()
ax.yaxis.offsetText.set_fontsize(fs)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.set_xlim(xmin=0, xmax=8.01)
#ax.set_ylim(ymin=0.993, ymax=1.001)
ax.set_ylim(ymin=0.835, ymax=1.01)
ax.legend(fontsize = fs)
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs-1)

ax = axs[1]

ax.semilogy(rho_dx05_t_new/ms_cu, np.abs(rho_dx05_max_new - rho_dx025_max_new), lw=lwd, linestyle='-', markersize=2, marker="x", markevery=mke, color='tab:blue', label=conv_lbl_dx05) #, alpha=0.7)
ax.semilogy(rho_dx025_t_new/ms_cu, 6.5*np.abs(rho_dx025_max_new - rho_dx0125_max_new), lw=lwd, linestyle='-', markersize=2, marker="x", markevery=mke, color='tab:green', label=conv_lbl_dx0125) #, alpha=0.7)

ax.set_xlabel('t [ms]', fontsize = fs)
ax.set_ylabel(r'$\epsilon$', fontsize = fs)
#ax = plt.axes()
ax.yaxis.offsetText.set_fontsize(fs)
plt.tick_params(axis='both', which='major', labelsize=fs)
ax.set_xlim(xmin=0, xmax=8.01)
ax.set_ylim(ymin=5e-6, ymax=5e-4)
#ax.set_ylim(ymin=0.85, ymax=1.01)
ax.legend(fontsize = fs)
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs-1)

fig.savefig(home+"/Desktop/magTOV_rhomax.pdf", bbox_inches = 'tight', pad_inches = 0.03)


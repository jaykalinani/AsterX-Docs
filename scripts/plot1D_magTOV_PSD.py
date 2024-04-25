#!/usr/bin/python

# Compute power spectral density of given data
#from plot_defaults import *
import matplotlib.pyplot as plt
from matplotlib.mlab import detrend_linear
import numpy as np
import matplotlib.ticker as mticker
import matplotlib as mpl
import os
home = os.environ["HOME"]

xlim = (0,12.1)
ylim_0p5 = (-133,-40)
ylim_0p25 = (-110,-47)
ylim_0p125 = (-147,-57)
t_start = 0.0
#t_end = 7.205

#mpl.rcParams['font.family'] = 'DejaVu Serif'
#mpl.rcParams['font.size'] = 24
#plt.rcParams['font.weight'] = 'bold'
#mpl.rcParams['axes.linewidth'] = 2
#mpl.rcParams['pdf.fonttype'] = 42
#mpl.rcParams['ps.fonttype'] = 42
#mpl.rcParams['mathtext.fontset'] = 'cm'
fs = 18
fontsize = 18
line_width = 3

# load data for GRAMX
#-----------------------------------------------------------------------------
Fx_0p5 = np.loadtxt("../TOV_data_GRaMX/paper_delta0p5_time.txt", unpack=True)
Fy_0p5 = np.loadtxt("../TOV_data_GRaMX/paper_delta0p5_rho.txt", unpack=True)
dt_0p5 = 1.0/(Fx_0p5[1]-Fx_0p5[0])

Fx_0p25 = np.loadtxt("../TOV_data_GRaMX/paper_delta0p25_time.txt", unpack=True)
Fy_0p25 = np.loadtxt("../TOV_data_GRaMX/paper_delta0p25_rho.txt", unpack=True)
dt_0p25 = 1.0/(Fx_0p25[1]-Fx_0p25[0])

Fx_0p125 = np.loadtxt("../TOV_data_GRaMX/paper_delta0p125_time.txt", unpack=True)
Fy_0p125 = np.loadtxt("../TOV_data_GRaMX/paper_delta0p125_rho.txt", unpack=True)
dt_0p125 = 1.0/(Fx_0p125[1]-Fx_0p125[0])
#-----------------------------------------------------------------------------


# load data for AsterX
#
dr_dx05 = "/Users/jaykalinani/Doc/PhD_Work/13_CarpetX/TOV_data_AsterX/magTOV_Z4c_AMR_dx05/"
dr_dx025 = "/Users/jaykalinani/Doc/PhD_Work/13_CarpetX/TOV_data_AsterX/magTOV_Z4c_AMR_dx025/"
dr_dx0125 = "/Users/jaykalinani/Doc/PhD_Work/13_CarpetX/TOV_data_AsterX/magTOV_Z4c_AMR_dx0125/"


rho_dx05 = np.genfromtxt(dr_dx05+"hydrobasex-rho.tsv")
rho_dx05_t = rho_dx05[:,1]/2.0301291587112965e+02
rho_dx05_max = rho_dx05[:,3]/rho_dx05[0,3]
dt_0p5_Ax = 1.0/(rho_dx05_t[1]-rho_dx05_t[0])

lbl_dx05 = r'$\Delta x_f=0.5$'

rho_dx025 = np.genfromtxt(dr_dx025+"hydrobasex-rho-all.tsv")
rho_dx025_t = rho_dx025[:,1]/2.0301291587112965e+02
rho_dx025_max = rho_dx025[:,3]/rho_dx025[0,3]
dt_0p25_Ax = 1.0/(rho_dx025_t[1]-rho_dx025_t[0])
lbl_dx025 = r'$\Delta x_f=0.25$'

rho_dx0125 = np.genfromtxt(dr_dx0125+"hydrobasex-rho.tsv")
rho_dx0125_t = rho_dx0125[:,1]/2.0301291587112965e+02
rho_dx0125_max = rho_dx0125[:,3]/rho_dx0125[0,3]
dt_0p125_Ax = 1.0/(rho_dx0125_t[1]-rho_dx0125_t[0])
lbl_dx0125 = r'$\Delta x_f=0.125$'
#

# mode names and frequencies
modes = {
"F" : 1.44252691528028,
"H1": 3.95408396916149,
"H2": 5.91494894170517,
"H3": 7.77382493405710,
"H4": 9.58768293806276,
"H5": 11.3772129983097,
#"H6": 13.1520241905666,
#"H7": 14.9172321735655,
}


# Plot basics
fig, ax = plt.subplots( 3, 1, figsize=(6.5,6.55), sharex='col', sharey='row', gridspec_kw={'hspace': 0.0, 'wspace': 0})

fig, ax = plt.subplots(3, 1,
            figsize=(8,8),
            sharex='col',
            sharey = 'row',
            #tight_layout=True,
            layout="constrained",
#            height_ratios = [1.8,1],
          #  width_ratios = [1,1],
          #  subplot_kw=dict(box_aspect=0.47),
            )

#fig.subplots_adjust(top=0.92,bottom=0.16, left=0.15,right=0.98)
#ax = fig.add_subplot(1,1,1)
lwd = 1.8
# plot modes as vertical lines and label them
for mode, freq in modes.items():
  ax[0].plot((freq,freq), ylim_0p5, color='k', linestyle='dotted', linewidth=lwd, alpha=0.6)
  ax[1].plot((freq,freq), ylim_0p25, color='k', linestyle='dotted', linewidth=lwd, alpha=0.6)
  ax[2].plot((freq,freq), ylim_0p125, color='k', linestyle='dotted', linewidth=lwd, alpha=0.6)
  ax[0].text(freq-0.1, ylim_0p5[1]+4, mode, fontsize=fontsize)

# plot PSD
#-----------------------------------------------------------------------------
frac_0p5 = 0.5 * len(Fy_0p5)
frac_0p5_Ax = 0.5 * len(rho_dx05_max)
#ax[0].psd(Fy_0p5, NFFT=int(frac_0p5), pad_to=5*int(frac_0p5), noverlap=int(frac_0p5*0.5), detrend=detrend_linear, Fs=dt_0p5, linestyle='solid', color='darkblue', scale_by_freq=False, linewidth=2.5)
ax[0].psd(rho_dx05_max, NFFT=int(frac_0p5_Ax), pad_to=5*int(frac_0p5_Ax), noverlap=int(frac_0p5_Ax*0.5), detrend=detrend_linear, Fs=dt_0p5_Ax, linestyle='solid', color='tab:blue', scale_by_freq=False, linewidth=lwd)

frac_0p25 = 0.5 * len(Fy_0p25)
frac_0p25_Ax = 0.5 * len(rho_dx025_max)
#ax[1].psd(Fy_0p25, NFFT=int(frac_0p25), pad_to=5*int(frac_0p25), noverlap=int(frac_0p25*0.5), detrend=detrend_linear, Fs=dt_0p25, linestyle='solid', color='darkgreen', scale_by_freq=False, linewidth=2.5)
ax[1].psd(rho_dx025_max, NFFT=int(frac_0p25_Ax), pad_to=5*int(frac_0p25_Ax), noverlap=int(frac_0p25_Ax*0.5), detrend=detrend_linear, Fs=dt_0p25_Ax, linestyle='solid', color='tab:green', scale_by_freq=False, linewidth=lwd)


frac_0p125 = 0.5 * len(Fy_0p125)
frac_0p125_Ax = 0.5 * len(rho_dx0125_max)
#ax[2].psd(Fy_0p125, NFFT=int(frac_0p125), pad_to=5*int(frac_0p125), noverlap=int(frac_0p125*0.5), detrend=detrend_linear, Fs=dt_0p125, linestyle='-', color='darkred', scale_by_freq=False, linewidth=2.5)
ax[2].psd(rho_dx0125_max, NFFT=int(frac_0p125_Ax), pad_to=5*int(frac_0p125_Ax), noverlap=int(frac_0p125_Ax*0.5), detrend=detrend_linear, Fs=dt_0p125_Ax, linestyle='solid', color='tab:red', scale_by_freq=False, linewidth=lwd)

#-----------------------------------------------------------------------------

# plot properties
#ax[0].set_xlim(xlim)

ax[0].set_xlim(xlim)
ax[0].set_ylim(ylim_0p5)
ax[1].set_xlim(xlim)
ax[1].set_ylim(ylim_0p25)
ax[2].set_xlim(xlim)
ax[2].set_ylim(ylim_0p125)

for i in range(3):
#    ax[i].tick_params(axis='y', labelsize=19 )
#    ax[i].tick_params(axis='x', labelsize=19 )
    ax[i].tick_params(axis='both', which='major', labelsize=fontsize-5)
#    ax[i].xaxis.set_tick_params(which='major', size=10, width=2, pad=8)
#    ax[i].xaxis.set_tick_params(which='minor', size=5, width=2)
#    ax[i].yaxis.set_tick_params(which='major', size=10, width=2, pad=8)
#    ax[i].yaxis.set_tick_params(which='minor', size=5, width=2, ) 

#ax[0].tick_params(axis='both', which='major', labelsize=fs-8)
#ax[1].tick_params(axis='both', which='major', labelsize=fs-8)
#ax[2].tick_params(axis='both', which='major', labelsize=fs-8)

#plt.title('PSD', fontsize=fontsize)
ax[0].set_xlabel('', fontsize=fontsize-3)
ax[0].set_ylabel('', fontsize=fontsize)
ax[1].set_xlabel('', fontsize=fontsize-3)
ax[1].set_ylabel(r'PSD [dB/Hz]', fontsize=fontsize)
ax[2].set_xlabel(r'$f$ [kHz]', fontsize=fontsize-3)
ax[2].set_ylabel('', fontsize=fontsize)
#fig.text(0.02, 0.5, r'PSD [dB/Hz]', va='center', rotation='vertical', fontsize=fontsize-3)

props = dict(boxstyle='round', facecolor='white', alpha=0.9)
ax[0].text(9.7, -57.5, r"$\Delta x_f=0.5$", fontsize=fontsize-6)
ax[1].text(9.7, -60, r"$\Delta x_f=0.25$", fontsize=fontsize-6)
ax[2].text(9.7, -75, r"$\Delta x_f=0.125$", fontsize=fontsize-6)

for i in range(3):
    ax[i].xaxis.set_major_locator(mpl.ticker.MultipleLocator(1.0))
    ax[i].xaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax[i].xaxis.grid(False)
    ax[i].yaxis.set_major_locator(mpl.ticker.MultipleLocator(25))
    ax[i].yaxis.set_minor_locator(mpl.ticker.MultipleLocator(12.5))
    ax[i].yaxis.grid(False)
    #set_tick_sizes(ax, 8, 4)

#plt.tight_layout()
#plt.savefig('all_PSD.pdf', format="pdf")
#plt.savefig('all_PSD.png')
#plt.show()
fig.savefig(home+"/Desktop/magTOV_PSD.pdf", bbox_inches = 'tight', pad_inches = 0.03)


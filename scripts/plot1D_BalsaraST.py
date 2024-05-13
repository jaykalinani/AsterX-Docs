#01_Balsara

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from numpy import *
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter

import os
home = os.environ["HOME"]

fname_B1 = "/Users/jaykalinani/Ubuntu_files/Balsara1_shocktube_xdir/Balsara1_shocktube_xdir/"
exact_data1 = np.genfromtxt("/Users/jaykalinani/Doc/PhD_Work/6_SpritzC2P/tests/scripts/ExactSolver_Giacomazzo/exact_solutions/Balsara1/solution_last.dat")

fname_B2 = "/Users/jaykalinani/Ubuntu_files/Balsara2_shocktube_xdir/Balsara2_shocktube_xdir/"
exact_data2 = np.genfromtxt("/Users/jaykalinani/Doc/PhD_Work/6_SpritzC2P/tests/scripts/ExactSolver_Giacomazzo/exact_solutions/Balsara2/solution_last.dat")

fname_B3 = "/Users/jaykalinani/Ubuntu_files/Balsara3_shocktube_xdir/Balsara3_shocktube_xdir/"
exact_data3 = np.genfromtxt("/Users/jaykalinani/Doc/PhD_Work/6_SpritzC2P/tests/scripts/ExactSolver_Giacomazzo/exact_solutions/Balsara3/solution_last.dat")

fname_B4 = "/Users/jaykalinani/Ubuntu_files/Balsara4_shocktube_xdir/Balsara4_shocktube_xdir/"
exact_data4 = np.genfromtxt("/Users/jaykalinani/Doc/PhD_Work/6_SpritzC2P/tests/scripts/ExactSolver_Giacomazzo/exact_solutions/Balsara4/solution_last.dat")

fname_B5 = "/Users/jaykalinani/Ubuntu_files/Balsara5_shocktube_xdir/Balsara5_shocktube_xdir/"
exact_data5 = np.genfromtxt("/Users/jaykalinani/Doc/PhD_Work/6_SpritzC2P/tests/scripts/ExactSolver_Giacomazzo/exact_solutions/Balsara5_psi/solution_last.dat")

tx1 = exact_data1[:,0]
rho1 = exact_data1[:,1]
by1 = exact_data1[:,6]

tx2 = exact_data2[:,0]
rho2 = exact_data2[:,1]
by2 = exact_data2[:,6]

tx3 = exact_data3[:,0]
rho3 = exact_data3[:,1]
by3 = exact_data3[:,6]

tx4 = exact_data4[:,0]
rho4 = exact_data4[:,1]
by4 = exact_data4[:,6]

tx5 = exact_data5[:,0]
rho5 = exact_data5[:,1]
by5 = exact_data5[:,6]

num = '002560'
num2 = '003520' #for Balsara 5

indi = -1936
indiB3 = -1582 
indj = -4 

lrpa='AsterX'

lwd=1
mkr1="x"
mkr2="+"
mkrs1 = 4
mkrs2 = 6
ml = MultipleLocator(5)

fs= 14.5
fig, axs = plt.subplots(5, 2, figsize=(13, 11.5), sharex=True) # constrained_layout=True)

rhox_B1 = np.genfromtxt(fname_B1+"hydrobasex-rho.it"+str(num)+".x.tsv")
by_B1 = np.genfromtxt(fname_B1+"hydrobasex-bvec.it"+str(num)+".x.tsv")

rhox_B2 = np.genfromtxt(fname_B2+"hydrobasex-rho.it"+str(num)+".x.tsv")
by_B2 = np.genfromtxt(fname_B2+"hydrobasex-bvec.it"+str(num)+".x.tsv")

rhox_B3 = np.genfromtxt(fname_B3+"hydrobasex-rho.it"+str(num)+".x.tsv")
by_B3 = np.genfromtxt(fname_B3+"hydrobasex-bvec.it"+str(num)+".x.tsv")

rhox_B4 = np.genfromtxt(fname_B4+"hydrobasex-rho.it"+str(num)+".x.tsv")
by_B4 = np.genfromtxt(fname_B4+"hydrobasex-bvec.it"+str(num)+".x.tsv")

rhox_B5 = np.genfromtxt(fname_B5+"hydrobasex-rho.it"+str(num2)+".x.tsv")
by_B5 = np.genfromtxt(fname_B5+"hydrobasex-bvec.it"+str(num2)+".x.tsv")

clr1 = 'r'
clr2 = 'g'

##B1

#rho
ax = axs[0,0]
ax.plot(rhox_B1[indi:indj,7], (rhox_B1[indi:indj,10]), lw=lwd, linestyle='None', markersize=mkrs1, marker=mkr1, color=clr1, label=lrpa)
ax.plot(tx1, (rho1), '-', color='k', label='Exact')
ax.set_ylabel(r'$\rho$', fontsize = fs)
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs)
ax.tick_params(direction='in', which='both', labelsize=fs-2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(plt.MaxNLocator(5))

ax.set_xlim(xmin=-0.5, xmax=0.5)
ax.legend(fontsize=fs-2)
#ax.set_ylim(ymin=-1.001, ymax=-0.994)

plt.text( -0.2, 0.5, r"$Balsara\,1$", color='black', rotation='vertical',fontsize=fs+3, transform = ax.transAxes,
          ha="center", va="center",bbox=dict(boxstyle="round",ec=(0.0, 0.0, 0.0),fc=(0.5, 0.5, 0.5), alpha=0.2) )
ax.grid(alpha=0.4)

#By
ax = axs[0,1]
ax.plot(by_B1[indi:indj,7], by_B1[indi:indj,11], lw=lwd, linestyle='None', markersize=mkrs1, marker=mkr1, color=clr1, label=lrpa)
ax.plot(tx1, by1, '-', color='k', label='Exact')
#ax.set_xlabel('x', fontsize = fs)
ax.set_ylabel(r'${\rm B^y}$', fontsize = fs)
ax.set_xlim(xmin=-0.5, xmax=0.5)
#ax.legend()
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs)
ax.tick_params(direction='in', which='both', labelsize=fs-2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
ax.grid(alpha=0.4)

####

##B2

#rho
ax = axs[1,0]
ax.plot(rhox_B2[indi:indj,7], (rhox_B2[indi:indj,10]), lw=lwd, linestyle='None', markersize=mkrs1, marker=mkr1, color=clr1, label=lrpa)
ax.plot(tx2, (rho2), '-', color='k', label='Exact')
ax.set_ylabel(r'$\rho$', fontsize = fs)
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs)
ax.tick_params(direction='in', which='both', labelsize=fs-2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
ax.grid(alpha=0.4)
ax.set_xlim(xmin=-0.5, xmax=0.5)

plt.text( -0.2, 0.5, r"$Balsara\,2$", color='black', rotation='vertical',fontsize=fs+3, transform = ax.transAxes,
          ha="center", va="center",bbox=dict(boxstyle="round",ec=(0.0, 0.0, 0.0),fc=(0.5, 0.5, 0.5), alpha=0.2) )

#By
ax = axs[1,1]
ax.plot(by_B2[indi:indj,7], by_B2[indi:indj,11], lw=lwd, linestyle='None', markersize=mkrs1, marker=mkr1, color=clr1, label=lrpa)
ax.plot(tx2, by2, '-', color='k', label='Exact')
ax.set_ylabel(r'${\rm B^y}$', fontsize = fs)
ax.set_xlim(xmin=-0.5, xmax=0.5)
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs)
ax.tick_params(direction='in', which='both', labelsize=fs-2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
ax.grid(alpha=0.4)


##B3

#rho
ax = axs[2,0]
ax.plot(rhox_B3[indiB3:indj,7], (rhox_B3[indiB3:indj,10]), lw=lwd, linestyle='None', markersize=mkrs1, marker=mkr1, color=clr1, label=lrpa)
ax.plot(tx3, (rho3), '-', color='k', label='Exact')
ax.set_ylabel(r'$\rho$', fontsize = fs)
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs)
ax.tick_params(direction='in', which='both', labelsize=fs-2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
ax.grid(alpha=0.4)
ax.set_xlim(xmin=-0.5, xmax=0.5)

plt.text( -0.2, 0.5, r"$Balsara\,3$", color='black', rotation='vertical',fontsize=fs+3, transform = ax.transAxes,
          ha="center", va="center",bbox=dict(boxstyle="round",ec=(0.0, 0.0, 0.0),fc=(0.5, 0.5, 0.5), alpha=0.2) )

#By
ax = axs[2,1]
ax.plot(by_B3[indiB3:indj,7], by_B3[indiB3:indj,11], lw=lwd, linestyle='None', markersize=mkrs1, marker=mkr1, color=clr1, label=lrpa)
ax.plot(tx3, by3, '-', color='k', label='Exact')
ax.set_ylabel(r'${\rm B^y}$', fontsize = fs)
ax.set_xlim(xmin=-0.5, xmax=0.5)
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs)
ax.tick_params(direction='in', which='both', labelsize=fs-2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
ax.grid(alpha=0.4)


##B4

#rho
ax = axs[3,0]
ax.plot(rhox_B4[indi:indj,7], (rhox_B4[indi:indj,10]), lw=lwd, linestyle='None', markersize=mkrs1, marker=mkr1, color=clr1, label=lrpa)
ax.plot(tx4, (rho4), '-', color='k', label='Exact')
ax.set_ylabel(r'$\rho$', fontsize = fs)
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs)
ax.tick_params(direction='in', which='both', labelsize=fs-2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(plt.MaxNLocator(5))

ax.set_xlim(xmin=-0.5, xmax=0.5)
ax.grid(alpha=0.4)

plt.text( -0.2, 0.5, r"$Balsara\,4$", color='black', rotation='vertical',fontsize=fs+3, transform = ax.transAxes,
          ha="center", va="center",bbox=dict(boxstyle="round",ec=(0.0, 0.0, 0.0),fc=(0.5, 0.5, 0.5), alpha=0.2) )
#By
ax = axs[3,1]
ax.plot(by_B4[indi:indj,7], by_B4[indi:indj,11], lw=lwd, linestyle='None', markersize=mkrs1, marker=mkr1, color=clr1, label=lrpa)
ax.plot(tx4, by4, '-', color='k', label='Exact')
ax.set_ylabel(r'${\rm B^y}$', fontsize = fs)
ax.set_xlim(xmin=-0.5, xmax=0.5)
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs)
ax.tick_params(direction='in', which='both', labelsize=fs-2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
ax.grid(alpha=0.4)


##B5

#rho
ax = axs[4,0]
ax.plot(rhox_B5[indi:indj,7], (rhox_B5[indi:indj,10]), lw=lwd, linestyle='None', markersize=mkrs1, marker=mkr1, color=clr1, label=lrpa)
ax.plot(tx5, (rho5), '-', color='k', label='Exact')
ax.set_ylabel(r'$\rho$', fontsize = fs)
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs)
ax.tick_params(direction='in', which='both', labelsize=fs-2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
ax.grid(alpha=0.4)

ax.set_xlim(xmin=-0.5, xmax=0.5)
ax.set_xlabel('x', fontsize = fs)

plt.text( -0.2, 0.5, r"$Balsara\,5$", color='black', rotation='vertical',fontsize=fs+3, transform = ax.transAxes,
          ha="center", va="center",bbox=dict(boxstyle="round",ec=(0.0, 0.0, 0.0),fc=(0.5, 0.5, 0.5), alpha=0.2) )

#By
ax = axs[4,1]
ax.plot(by_B5[indi:indj,7], by_B5[indi:indj,11], lw=lwd, linestyle='None', markersize=mkrs1, marker=mkr1, color=clr1, label=lrpa)
ax.plot(tx5, by5, '-', color='k', label='Exact')
ax.set_ylabel(r'${\rm B^y}$', fontsize = fs)
ax.set_xlim(xmin=-0.5, xmax=0.5)
for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(fs)
ax.tick_params(direction='in', which='both', labelsize=fs-2)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_major_locator(plt.MaxNLocator(5))
ax.set_xlabel('x', fontsize = fs)
ax.grid(alpha=0.4)

fig.savefig(home+"/Desktop/1D_BalsaraTests.pdf", bbox_inches = 'tight', pad_inches = 0.1)


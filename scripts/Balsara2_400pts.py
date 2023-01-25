#01_Balsara

import numpy as np
from postcactus.simdir import SimDir
import matplotlib.pyplot as plt
from postcactus import grid_data as gd

import os
home = os.environ["HOME"]


fname_LR = "/Users/jaykalinani/Ubuntu_files/Balsara2_shocktube/"
fname_spr = "/Users/jaykalinani/Ubuntu_files/Balsara2_Spritz_Noble_LxF_PPM_400pts/"
#exact_data = np.genfromtxt("/Users/jaykalinani/Doc/PhD_Work/13_CarpetX/scripts/ExactSolver_Giacomazzo/rmhd_riemann_solver/solution.dat")
exact_data = np.genfromtxt("/Users/jaykalinani/Doc/PhD_Work/13_CarpetX/scripts/ExactSolver_Giacomazzo/exact_solutions/Balsara2/solution_last.dat")
tx = exact_data[:,0]
rhox = exact_data[:,1]
pt = exact_data[:,2]
by = exact_data[:,6]
bz = exact_data[:,7]
vx = exact_data[:,3]
vy = exact_data[:,4]
vz = exact_data[:,5]

num = '000640' #001664
num2 = '000640' #000800
indi = -404 #-22# -140 #-240
indj =  -4 #-11 #-4 #-4

xl = -.5
xr = 0.5 

lspr='Spritz: 400pt'
#lgrx='AsterX: 300pt'
lgrx='AsterX: 400pt'

lwd=1
mkr1="o"
mkr2="x"
mkrs = 3
rhox_lr = np.genfromtxt(fname_LR+"hydrobase-rho.it"+str(num2)+".x.tsv")
rhox_spr = np.genfromtxt(fname_spr+"rho.x.asc")

vel_lr = np.genfromtxt(fname_LR+"hydrobase-vel.it"+str(num2)+".x.tsv")
velx_spr = np.genfromtxt(fname_spr+"vel[0].x.asc")
vely_spr = np.genfromtxt(fname_spr+"vel[1].x.asc")
velz_spr = np.genfromtxt(fname_spr+"vel[2].x.asc")

eps_lr = np.genfromtxt(fname_LR+"hydrobase-eps.it"+str(num2)+".x.tsv")
eps_spr = np.genfromtxt(fname_spr+"eps.x.asc")

press_lr = np.genfromtxt(fname_LR+"hydrobase-press.it"+str(num2)+".x.tsv")
press_spr = np.genfromtxt(fname_spr+"press.x.asc")

#####rho_x plot
fig = plt.figure(tight_layout=True)
plt.plot(rhox_spr[indi:indj,9], rhox_spr[indi:indj,12], lw=lwd, linestyle='--', markersize=mkrs, marker=mkr1, color='magenta', label=lspr)
plt.plot(rhox_lr[:,7], rhox_lr[:,10], lw=lwd,linestyle='none',markersize=mkrs, marker=mkr2, color='r', label=lgrx)
plt.plot(tx, rhox, '-', color='k', label='Exact', lw=lwd)
plt.xlabel('x', fontsize = 16)
plt.ylabel(r'$\rho$', fontsize = 16)
ax = plt.axes()
ax.yaxis.offsetText.set_fontsize(16)
plt.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlim(xmin=xl, xmax=xr)
ax.legend()
fig.savefig(home+"/Desktop/rho_x.pdf")

##velx
fig = plt.figure(tight_layout=True)
plt.plot(velx_spr[indi:indj,9], velx_spr[indi:indj,12], lw=lwd,linestyle='--',markersize=mkrs, marker=mkr1, color='magenta', label=lspr)
plt.plot(vel_lr[:,7], vel_lr[:,10],lw=lwd,linestyle='-.',markersize=mkrs, marker=mkr2, color='r', label=lgrx)
plt.plot(tx, vx, '-', color='k',  label='Exact', lw=lwd)
plt.xlabel('x', fontsize = 16)
plt.ylabel(r'$\rm v_x$', fontsize = 16)
ax = plt.axes()
ax.yaxis.offsetText.set_fontsize(16)
plt.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlim(xmin=xl, xmax=xr)
ax.legend()
fig.savefig(home+"/Desktop/velx_x.pdf")

##vely
fig = plt.figure(tight_layout=True)
plt.plot(vely_spr[indi:indj,9], vely_spr[indi:indj,12], lw=lwd,linestyle='--',markersize=mkrs, marker=mkr1, color='magenta', label=lspr)
plt.plot(vel_lr[:,7], vel_lr[:,11], lw=lwd,linestyle='-.',markersize=mkrs, marker=mkr2, color='r', label=lgrx)
plt.plot(tx, vy, '-', color='k', label='Exact')
plt.xlabel('x', fontsize = 16)
plt.ylabel(r'$\rm v_y$', fontsize = 16)
ax = plt.axes()
ax.yaxis.offsetText.set_fontsize(16)
plt.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlim(xmin=xl, xmax=xr)
ax.legend()
fig.savefig(home+"/Desktop/vely_x.pdf")

##velz
fig = plt.figure(tight_layout=True)
plt.plot(velz_spr[indi:indj,9], velz_spr[indi:indj,12], lw=lwd,linestyle='--',markersize=mkrs, marker=mkr1, color='magenta', label=lspr)
plt.plot(vel_lr[:,7], vel_lr[:,12], lw=lwd,linestyle='-.',markersize=mkrs, marker=mkr2, color='r', label=lgrx)
plt.plot(tx, vz, '-', color='k', label='Exact')
plt.xlabel('x', fontsize = 16)
plt.ylabel(r'$\rm v_z$', fontsize = 16)
ax = plt.axes()
ax.yaxis.offsetText.set_fontsize(16)
plt.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlim(xmin=xl, xmax=xr)
ax.legend()
fig.savefig(home+"/Desktop/velz_x.pdf")

##eps
fig = plt.figure(tight_layout=True)
bx =np.empty(np.size(exact_data[:,1]))
bx.fill(5.0)
#bx.fill(0.0)
gamma =5.0/3.0
wx_comp = 1.0/np.sqrt(1.0 - (vx*vx+vy*vy+vz*vz))

term1 = (bx*bx + by*by + bz*bz)/(wx_comp*wx_comp)
term2 = (bx*vx + by*vy + bz*vz)**2
pm = (term1 + term2)/2.0
pg = pt - pm

eps_exact = pg/((gamma-1.0)*rhox)

plt.plot(eps_spr[indi:indj,9], eps_spr[indi:indj,12], lw=lwd,linestyle='--',markersize=mkrs, marker=mkr1, color='magenta', label=lspr)
plt.plot(eps_lr[:,7], eps_lr[:,10], lw=lwd,linestyle='-.',markersize=mkrs, marker=mkr2, color='r', label=lgrx)
plt.plot(exact_data[:,0], eps_exact, '-', color='k', label='Exact')
plt.xlabel('x', fontsize = 16)
plt.ylabel(r'$\epsilon$', fontsize = 16)
ax = plt.axes()
ax.yaxis.offsetText.set_fontsize(16)
ax.set_xlim(xmin=xl, xmax=xr)
plt.tick_params(axis='both', which='major', labelsize=16)
ax.legend()
fig.savefig(home+"/Desktop/eps_x.pdf")

#####press gas
fig = plt.figure(tight_layout=True)
plt.plot(press_spr[indi:indj,9], press_spr[indi:indj,12], lw=lwd,linestyle='--',markersize=mkrs, marker=mkr1, color='magenta', label=lspr)
plt.plot(press_lr[:,7], press_lr[:,10], lw=lwd,linestyle='none',markersize=mkrs, marker=mkr2, color='r', label=lgrx)
plt.plot(exact_data[:,0], pg, '-', color='k', label='Exact', lw=lwd)
plt.xlabel('x', fontsize = 16)
plt.ylabel(r'${\rm P_g}$', fontsize = 16)
ax = plt.axes()
ax.yaxis.offsetText.set_fontsize(16)
ax.set_xlim(xmin=xl, xmax=xr)
plt.tick_params(axis='both', which='major', labelsize=16)
ax.legend()
fig.savefig(home+"/Desktop/Pg_x.pdf")



#import matplotlib.backends.backend_pdf
#pdf = matplotlib.backends.backend_pdf.PdfPages(home+"/Desktop/Balsara1.pdf")
#for figx in xrange(1, fig().number): ## will open an empty extra figure :(
#    pdf.savefig( figx )
#pdf.close()



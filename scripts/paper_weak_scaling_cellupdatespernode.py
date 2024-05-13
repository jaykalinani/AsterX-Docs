import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

SMALL_SIZE = 9
MEDIUM_SIZE = 9
BIGGER_SIZE = 9

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)
plt.rc('text', usetex=True)

fig = plt.figure(figsize=(6.5,4.5))
ax    = fig.add_subplot(111)
#plt.figure(figsize=(8.5,5.75))
ax.set_xscale('log') #, basex=10)
ax.set_yscale('log') #, basey=10)

nodes = np.array([8, 27, 64, 216, 512, 1728, 4096])
cores = nodes*8
zcpersec_CowUni128 = np.array(
[
8.55,
25.3,
57.1,
178,
413,
1230,
2940
] ) 

zcpersec_CowUni240 = np.array(
[
10,
31.3,
66.7,
229,
530,
1700,
3980
] )


zcpersec_Z4cUni128 = np.array(
[
4.1,
12.6,
29.4,
94.3,
213,
439,
490
] )

zcpersec_Z4cUni240 = np.array(
[
4.62,
14.7,
33.5,
114,
261,
711,
1580
] )

zcpersec_Z4cAMR128= np.array(
[
0.954,
2.79,
5.90,
18.7,
24.4,
32.6,
170
] )

zcpersec_Z4cAMR240 = np.array(
[
1.56,
4.22,
7.95,
26,
43.9,
137,
179
] )

zcpersecpergpu_CowUni128 = zcpersec_CowUni128/cores
zcpersecpergpu_CowUni240 = zcpersec_CowUni240/cores

zcpersecpergpu_Z4cUni128 = zcpersec_Z4cUni128/cores
zcpersecpergpu_Z4cUni240 = zcpersec_Z4cUni240/cores

zcpersecpergpu_Z4cAMR128 = zcpersec_Z4cAMR128/cores
zcpersecpergpu_Z4cAMR240 = zcpersec_Z4cAMR240/cores
#idealrate = (cores/cores[0]) * zcpersec [0]
#ax.plot(cores, idealrate, '-^', color='darkblue', linewidth=1.1, label='Ideal')

ax.plot(cores, zcpersecpergpu_CowUni128
,'-o', 
color='red', linewidth=1.1, 
label=r'$SetupA: 128^3$ cells/GPU', alpha=0.7)

ax.plot(cores, zcpersecpergpu_Z4cUni128
,'-o',
color='green', linewidth=1.1,
label=r'$SetupB: 128^3$ cells/GPU', alpha=0.7)

ax.plot(cores, zcpersecpergpu_Z4cAMR128
,'-o',
color='blue', linewidth=1.1,
label=r'$SetupC: 128^3$ cells/GPU', alpha=0.7)

ax.plot(cores, zcpersecpergpu_CowUni240
,'--o',
color='darkred', linewidth=1.1,
label=r'$SetupA: 240^3$ cells/GPU', alpha=0.9)

ax.plot(cores, zcpersecpergpu_Z4cUni240
,'--o',
color='darkgreen', linewidth=1.1,
label=r'$SetupB: 240^3$ cells/GPU', alpha=0.9)

ax.plot(cores, zcpersecpergpu_Z4cAMR240
,'--o',
color='darkblue', linewidth=1.1,
label=r'$SetupC: 240^3$ cells/GPU', alpha=0.9)


ax.set_xlabel(r'Number of OLCF Frontier GPUs', fontsize=14)
ax.set_ylabel(r'zone-cycles / sec / GPU ($\times 10^8$) ', fontsize=14)
ax.tick_params(direction='in', which='both')
#ax.set_title("Magnetized TOV (AsterX + Cowling + idealgas + uniform grid)", fontsize=14)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,1))
#.grid(alpha=0.4)

plt.text(50, 0.3, r'Weak Scaling', fontsize=16)
ax.legend(fontsize=10.5, ncol= 2, loc='lower left') 

for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(13)
#plt.ylim(8e6, 2e9)


ax.set_xlim(xmin=40, xmax=50000) #xmax=155)   #xmax=data[1][1].t[-1] / MS_CU)
ax.set_ylim(ymin=5e-4, ymax=5e-1)


plt.savefig("/Users/jaykalinani/Desktop/weak_scaling.pdf", bbox_inches = 'tight', pad_inches = 0.03)
#plt.savefig("/Users/jaykalinani/Desktop/cellupdates_strong_scaling.png", bbox_inches = 'tight', pad_inches = 0.03)

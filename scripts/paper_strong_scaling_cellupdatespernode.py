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

#fig = plt.figure(figsize=(6.5,3.75))
fig = plt.figure(figsize=(5.5,4.5))
ax    = fig.add_subplot(111)
#plt.figure(figsize=(8.5,5.75))
ax.set_xscale('log') #, basex=10)
ax.set_yscale('log') #, basey=10)
#plt.xticks(np.array([8, 16, 32, 64, 128, 256, 512, 1024]), ["8", "16", "32", "64", "128", "256", "512", "1024"])

##########OLD######
#plt.plot(cores, idealrate, '-^', color='darkblue', linewidth=1.1)
#plt.plot(np.array([8, 16, 32, 64, 128]), np.array([0.0240679, 0.0240679*2,  0.0240679*4, 0.0240679*8, 0.0240679*16 ]),'-^', color='blue', linewidth=1.1, label=r'Ideal Rate')

#plt.plot(np.array([8, 16, 32, 64, 128, 256]), np.array([16049644.4253103, 26749861.2990325, 40223728.6818827, 58296505.854607, 74520693.3721593,  87153244.7937287]),'-o', color='blue', linewidth=1.1, label=r'IGM + ML_BSSN + Subcycling (Frontera)', alpha=0.7)

#plt.plot(np.array([8, 16, 32, 64, 128, 256, 512]), np.array([3.21851e+06, 4.86786e+06, 9.03244e+06, 1.25005e+07, 1.84129e+07,  2.12953e+07, 2.4786e+07]),'-o', color='green', linewidth=1.1, label=r'AsterX + Z4c + No Subcycling (Frontera)', alpha=0.7)

#plt.plot(np.array([8, 16, 32]), np.array([3.5326e+06, 6.20103e+06, 9.29504e+06 ]),'-o', color='orange', linewidth=1.1, label=r'AsterX + Z4c + No Subcycling (Frontera, Optimized)', alpha=0.7)

#plt.plot(np.array([8, 16, 32, 64, 128, 256, 512]), np.array([7.32585e+06, 1.35079e+07, 2.33422e+07, 3.68688e+07, 5.23376e+07, 6.00127e+07, 7.8675e+07]),'-o', color='red', linewidth=1.1, label=r'AsterX + Z4c + No Subcycling (Frontier)', alpha=0.7)

#plt.plot(np.array([8, 16, 32, 64, 128, 256, 512, 1024]), np.array([ 1.7638e+07, 3.04553e+07, 4.81173e+07, 7.08822e+07, 9.01771e+07, 1.02087e+08, 1.20474e+08, 1.52415e+08]),'-o', color='orange', linewidth=1.1, label=r'AsterX + Minkowski + No Subcycling (Frontier)', alpha=0.7)

#plt.ticklabel_format(axis='y', style='sci', scilimits=None, useOffset=None, useLocale=None, useMathText=True)
#######



nodes = np.array([8, 27, 64, 125, 216, 512, 1000])
cores = nodes*8
zcpersec_CowUni = np.array(
[
10.2,
26.9,
43.5,
66.5,
96.8,
150,
165
])

zcpersec_Z4cUni = np.array(
[
4.62,
13.5,
24.1,
37.3,
49.4,
76.6,
115
])

zcpersec_Z4cAMR = np.array(
[
1.56,
3.53,
5.50,
7.45,
8.62,
11.4,
10.9
])


#idealrate = (cores/cores[0]) * zcpersec [0]
#ax.plot(cores, idealrate, '-^', color='darkblue', linewidth=1.1, label='Ideal')

#128
ax.plot(cores, zcpersec_CowUni 
,'-o', 
color='red', linewidth=1.1, 
label=r'$SetupA$', alpha=0.7)

ax.plot(cores, zcpersec_Z4cUni
,'-o',
color='green', linewidth=1.1,
label=r'$SetupB$', alpha=0.7)

ax.plot(cores, zcpersec_Z4cAMR
,'-o',
color='blue', linewidth=1.1,
label=r'$SetupC$', alpha=0.7)

ax.set_xlabel(r'Number of OLCF Frontier GPUs', fontsize=14)
ax.set_ylabel(r'zone-cycles / sec ($\times 10^8$) ', fontsize=14)
ax.tick_params(direction='in', which='both')
plt.text(55, 200, r'Strong Scaling', fontsize=16)


#ax.set_title("Magnetized TOV (AsterX + Cowling + idealgas + uniform grid)", fontsize=14)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(-1,1))
#.grid(alpha=0.4)
ax.legend(fontsize=11.5, loc='lower right') #bbox_to_anchor = (0.6, 0.6), loc='center right') 

for item in (ax.get_yticklabels() + ax.get_xticklabels()):
      item.set_fontsize(13)
#plt.ylim(8e6, 2e9)


ax.set_xlim(xmin=45, xmax=11000) #xmax=155)   #xmax=data[1][1].t[-1] / MS_CU)
ax.set_ylim(ymin=1, ymax=300)


plt.savefig("/Users/jaykalinani/Desktop/strong_scaling.pdf", bbox_inches = 'tight', pad_inches = 0.03)
#plt.savefig("/Users/jaykalinani/Desktop/cellupdates_strong_scaling.png", bbox_inches = 'tight', pad_inches = 0.03)

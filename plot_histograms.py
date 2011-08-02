#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import string

mpl.rcParams['text.usetex'    ] = True
mpl.rcParams['font.size'      ] = 14
mpl.rcParams['axes.labelsize' ] = 14
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14 
mpl.rcParams['text.latex.unicode']=True

lpath=sys.argv[1]
spath=sys.argv[2]

#lpath='/home/james/code/auswertung_galaxy_files/data/'
#spath='/home/james/code/auswertung_galaxy_files/'

# Determine which histograms should be plotted:
histograms=('stf','mcold','rdisk','coldmetal','velocity', \
            'velocity_x','velocity_y','velocity_z',       \
            'pos_x', 'pos_y', 'pos_z', 'distance')
#histograms=('mcold')        

# Determine which other plots should be made:
#plots=('pos3d')
plots=('pos3d','cluster_stf','cluster_stf_avg','coldmetal_avg')
#plots=('')

# Get number of timesteps
tsteps=np.loadtxt(lpath+'/info.dat')[0]
ngalxs=np.loadtxt(lpath+'/info.dat')[1]

# Determine which timestep should be plotted
tstart=1286
#tend=sum(1 for line in open(lpath+'coldmetal_avg.dat')) #number of galaxy files
tend=1286


if 'stf' in histograms:
    fname=lpath+'/stf_histo.dat'
    histo=np.fromfile(fname,dtype='float32')
    histo=histo.reshape(tsteps,ngalxs+2)
    histo=np.hsplit(histo,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(histo[1][tstep-1][:], 40, range=(0,60),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(1,100000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ star\ formation}$')
        ax1.set_xlabel(r'$\mathrm{star\ formation\ (}\mathrm{m}_{\odot}/\mathrm{year)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%histo[0][tstep-1][1])
        fig.savefig(spath+'stf_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'stf_histo_%04i.pdf'%int(i))



if 'mcold' in histograms:
    fname=lpath+'/mcold_histo.dat'
    histo=np.fromfile(fname,dtype='float32')
    histo=histo.reshape(tsteps,ngalxs+2)
    histo=np.hsplit(histo,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(histo[1][tstep-1][:], 40, range=(0,50),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(1,100000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ gas\ disk\ mass}$')
        ax1.set_xlabel(r'$\mathrm{gas\ disk\ mass\ (}10^{10}\mathrm{m}_{\odot}/\mathrm{year)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%histo[0][tstep-1][1])
        fig.savefig(spath+'mcold_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'mcold_histo_%04i.pdf'%int(i))


if 'rdisk' in histograms:
    fname=lpath+'/rdisk_histo.dat'
    histo=np.fromfile(fname,dtype='float32')
    histo=histo.reshape(tsteps,ngalxs+2)
    histo=np.hsplit(histo,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(histo[1][tstep-1][:], 40, range=(0,50),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(1,100000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ disk\ scale\ length}$')
        ax1.set_xlabel(r'$\mathrm{disk\ scale\ length\ (kpc)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%histo[0][tstep-1][1])
        fig.savefig(spath+'rdisk_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'rdisk_histo_%04i.pdf'%int(i))


if 'coldmetal' in histograms:
    fname=lpath+'/coldmetal_histo.dat'
    histo=np.fromfile(fname,dtype='float32')
    histo=histo.reshape(tsteps,ngalxs+2)
    histo=np.hsplit(histo,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(histo[1][tstep-1][:], 40, range=(0,3),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(1,10000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ disk\ gas\ metallicity}$')
        ax1.set_xlabel(r'$\mathrm{gas\ disk\ metallicity\ (solar\ units)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%histo[0][tstep-1][1])
        fig.savefig(spath+'diskmetal_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'diskmetal_histo_%04i.pdf'%int(i))

# Attention, is velocity actually in comoving units? check!

if 'velocity' in histograms:
    fname=lpath+'/velocity_x.dat'
    vel_x=np.fromfile(fname,dtype='float32')
    vel_x=vel_x.reshape(tsteps,ngalxs+2)
    vel_x=np.hsplit(vel_x,[2])
    fname=lpath+'/velocity_y.dat'
    vel_y=np.fromfile(fname,dtype='float32')
    vel_y=vel_y.reshape(tsteps,ngalxs+2)
    vel_y=np.hsplit(vel_y,[2])
    fname=lpath+'/velocity_z.dat'
    vel_z=np.fromfile(fname,dtype='float32')
    vel_z=vel_z.reshape(tsteps,ngalxs+2)
    vel_z=np.hsplit(vel_z,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(np.sqrt(vel_x[1][tstep-1][:]**2+\
                                            vel_y[1][tstep-1][:]**2+\
                                            vel_z[1][tstep-1][:]**2), 40, range=(0,10000),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(1,10000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ total\ velocity}$')
        ax1.set_xlabel(r'$\mathrm{velocity\ (km/s)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%vel_x[0][tstep-1][1])
        fig.savefig(spath+'vel_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'vel_histo_%04i.pdf'%int(i))

if 'velocity_x' in histograms:
    fname=lpath+'/velocity_x.dat'
    histo=np.fromfile(fname,dtype='float32')
    histo=histo.reshape(tsteps,ngalxs+2)
    histo=np.hsplit(histo,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(histo[1][tstep-1][:], 40, range=(-5000,5000),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(1,10000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ x-velocity}$')
        ax1.set_xlabel(r'$\mathrm{velocity\ (km/s)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%histo[0][tstep-1][1])
        fig.savefig(spath+'vel_x_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'vel_x_histo_%04i.pdf'%int(i))
       
if 'velocity_y' in histograms:
    fname=lpath+'/velocity_y.dat'
    histo=np.fromfile(fname,dtype='float32')
    histo=histo.reshape(tsteps,ngalxs+2)
    histo=np.hsplit(histo,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(histo[1][tstep-1][:], 40, range=(-5000,5000),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(1,10000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ y-velocity}$')
        ax1.set_xlabel(r'$\mathrm{velocity\ (km/s)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%histo[0][tstep-1][1])
        fig.savefig(spath+'vel_y_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'vel_y_histo_%04i.pdf'%int(i))

if 'velocity_z' in histograms:
    fname=lpath+'/velocity_z.dat'
    histo=np.fromfile(fname,dtype='float32')
    histo=histo.reshape(tsteps,ngalxs+2)
    histo=np.hsplit(histo,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(histo[1][tstep-1][:], 40, range=(-5000,5000),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(1,10000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ z-velocity}$')
        ax1.set_xlabel(r'$\mathrm{velocity\ (km/s)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%histo[0][tstep-1][1])
        fig.savefig(spath+'vel_z_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'vel_z_histo_%04i.pdf'%int(i))

if 'pos_x' in histograms:
    fname=lpath+'/positions_x.dat'
    histo=np.fromfile(fname,dtype='float32')
    histo=histo.reshape(tsteps,ngalxs+2)
    histo=np.hsplit(histo,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(histo[1][tstep-1][:], 40, range=(-25000,25000),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(100,10000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ galaxies\ x-position}$')
        ax1.set_xlabel(r'$\mathrm{x-position\ (kpc)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%histo[0][tstep-1][1])
        fig.savefig(spath+'pos_x_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'pos_x_histo_%04i.pdf'%int(i))


if 'pos_y' in histograms:
    fname=lpath+'/positions_y.dat'
    histo=np.fromfile(fname,dtype='float32')
    histo=histo.reshape(tsteps,ngalxs+2)
    histo=np.hsplit(histo,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(histo[1][tstep-1][:], 40, range=(-25000,25000),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(100,10000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ galaxies\ y-position}$')
        ax1.set_xlabel(r'$\mathrm{y-position\ (kpc)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%histo[0][tstep-1][1])
        fig.savefig(spath+'pos_y_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'pos_y_histo_%04i.pdf'%int(i))

if 'pos_z' in histograms:
    fname=lpath+'/positions_z.dat'
    histo=np.fromfile(fname,dtype='float32')
    histo=histo.reshape(tsteps,ngalxs+2)
    histo=np.hsplit(histo,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(histo[1][tstep-1][:], 40, range=(-25000,25000),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(100,10000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ galaxies\ z-position}$')
        ax1.set_xlabel(r'$\mathrm{z-position\ (kpc)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%histo[0][tstep-1][1])
        fig.savefig(spath+'pos_z_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'pos_z_histo_%04i.pdf'%int(i))

if 'distance' in histograms:
    fname=lpath+'/positions_x.dat'
    pos_x=np.fromfile(fname,dtype='float32')
    pos_x=pos_x.reshape(tsteps,ngalxs+2)
    pos_x=np.hsplit(pos_x,[2])
    fname=lpath+'/positions_y.dat'
    pos_y=np.fromfile(fname,dtype='float32')
    pos_y=pos_y.reshape(tsteps,ngalxs+2)
    pos_y=np.hsplit(pos_y,[2])
    fname=lpath+'/positions_z.dat'
    pos_z=np.fromfile(fname,dtype='float32')
    pos_z=pos_z.reshape(tsteps,ngalxs+2)
    pos_z=np.hsplit(pos_z,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        n, bins, patches = ax1.hist(np.sqrt(pos_x[1][tstep-1][:]**2+\
                                            pos_y[1][tstep-1][:]**2+\
                                            pos_z[1][tstep-1][:]**2), 40, range=(0,30000),\
                    log=True, normed=False, align='mid',histtype='stepfilled',\
                    facecolor='0.85',color='black')
        ax1.grid(True)
        ax1.set_ylim(1,10000)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ galaxy\ distance\ from\ origin}$')
        ax1.set_xlabel(r'$\mathrm{distance\ (kpc)}$')
        ax1.set_ylabel(r'$\mathrm{\#\ of\ galaxies}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%pos_x[0][tstep-1][1])
        fig.savefig(spath+'distance_histo_%04i.png'%int(tstep),dpi=120)
        # fig.savefig(spath+'distance_histo_%04i.pdf'%int(i))

if 'pos3d' in plots:
    fname=lpath+'/positions_x.dat'
    pos_x=np.fromfile(fname,dtype='float32')
    pos_x=pos_x.reshape(tsteps,ngalxs+2)
    pos_x=np.hsplit(pos_x,[2])
    fname=lpath+'/positions_y.dat'
    pos_y=np.fromfile(fname,dtype='float32')
    pos_y=pos_y.reshape(tsteps,ngalxs+2)
    pos_y=np.hsplit(pos_y,[2])
    fname=lpath+'/positions_z.dat'
    pos_z=np.fromfile(fname,dtype='float32')
    pos_z=pos_z.reshape(tsteps,ngalxs+2)
    pos_z=np.hsplit(pos_z,[2])
    for tstep in range(tstart,tend+1,1):
        fig = plt.figure()
        ax1 = Axes3D(fig)
        plimit=5000
        mask=np.zeros(int(ngalxs))
        mask=mask.astype(bool)
        for i in range(0,int(ngalxs)):
            mask[i]=(abs(pos_x[1][tstep-1][i])<plimit and \
                         abs(pos_y[1][tstep-1][i])<plimit and \
                         abs(pos_z[1][tstep-1][i])<plimit)
        ax1.scatter(pos_x[1][tstep-1][mask],pos_y[1][tstep-1][mask],pos_z[1][tstep-1][mask],alpha=0.75)
        ax1.set_xlim3d(-plimit,plimit)	 
        ax1.set_ylim3d(-plimit,plimit)
        ax1.set_zlim3d(-plimit,plimit)
        ax1.set_title(r'$\mathrm{Histogram\ of\ the\ galaxy\ distance\ from\ center}$')
        ax1.set_xlabel(r'$\mathrm{(kpc)}$')
        ax1.set_ylabel(r'$\mathrm{(kpc)}$')
        ax1.set_zlabel(r'$\mathrm{(kpc)}$')
        fig.text(0.80,0.95,r'$\mathrm{z = }$')
        fig.text(0.85,0.95,r'$\mathrm{%.2f}$'%pos_x[0][tstep-1][1])
        fig.savefig(spath+'pos3d_%04i.png'%int(tstep),dpi=120)
        #fig.savefig(spath+'pos3d_%04i.pdf'%int(tstep))

if 'cluster_stf' in plots:
    data=np.fromfile(lpath+'/cluster_stf.dat',dtype='float32')
    data=data.reshape(tsteps,3)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(data[:,0]-13.7,data[:,2],'k')
    ax1.set_xlabel(r'lookback-time')
    ax1.set_title(r'global\ star\ formation')
    ax1.set_ylabel(r'star\ formation\ [$M_\odot$/year?]')
    ax1.set_yscale('log')
    fig.savefig(spath+'cluster_stf.png')

if 'cluster_stf_avg' in plots:
    data=np.fromfile(lpath+'/cluster_stf_avg.dat',dtype='float32')
    data=data.reshape(tsteps,3)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(data[:,0]-13.7,data[:,2],'k')
    ax1.set_xlabel(r'lookback-time')
    ax1.set_ylabel(r'average\ star\ formation\ [$M_\odot$/year?]')
    ax1.set_title(r'average\ star\ formation')
    ax1.set_yscale('log')
    fig.savefig(spath+'cluster_stf_avg.png')


if 'coldmetal_avg' in plots:
    data=np.fromfile(lpath+'/coldmetal_avg.dat',dtype='float32')
    data=data.reshape(tsteps,3)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(data[:,0]-13.7,data[:,2],'k')
    ax1.set_xlabel(r'lookback-time')
    ax1.set_title(r'average\ cold\ metal?')
    fig.savefig(spath+'coldmetal_avg.png')

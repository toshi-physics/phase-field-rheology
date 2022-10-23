#!/bin/env python3
# plot_field_velocity.py

import sys

# This code must be run with python3!
if (sys.version_info < (3, 5)):
    print("This code must be run with Python version 3.5 or higher")
    sys.exit(1)

import numpy as np
import matplotlib as mpl
import matplotlib.cm as mplcm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc
from matplotlib.animation import FuncAnimation

args = sys.argv
if (len(args) != 23):
    print("Usage: plot_nematic_field.py npoints lx ly Dr dt xxcol yycol xycol val_min val_max tic_start tic_end tic_inc tstart tend tinc make_movie print_to_screen field_file_dir field_file_name pos_file out_file")
    sys.exit(1)

npoints = int(args.pop(1))
lx = int(args.pop(1))
ly = int(args.pop(1))
Dr = float(args.pop(1))
dt = float(args.pop(1))
xxcol = int(args.pop(1))
yycol = int(args.pop(1))
xycol = int(args.pop(1))
val_min = float(args.pop(1))
val_max = float(args.pop(1))
tic_start = float(args.pop(1))
tic_end = float(args.pop(1))
tic_inc = float(args.pop(1))
tstart = int(args.pop(1))
tend = int(args.pop(1))
tinc = int(args.pop(1))
make_movie = bool(int(args.pop(1)))
print_to_screen = bool(int(args.pop(1)))
field_file_dir = args.pop(1)
field_file_name = args.pop(1)
pos_file = args.pop(1)
out_file = args.pop(1)
nframes = (tend-tstart)//tinc+1

use_label = 1 # 1
use_cbar = 1 # 1
sign = 1 # 1 or -1
head = 0 # 1

if (not make_movie):
    tend = tstart

# Data arrays
vec_len = (0.5*np.sqrt(lx*ly/npoints))
pos = np.zeros((nframes,npoints,2))
#vec = np.zeros((lx,ly,3))
field = np.zeros((lx,ly))
xa = np.arange(lx)
ya = np.arange(ly)
xx, yy = np.meshgrid(xa,ya)
vx = np.zeros((lx,ly))
vy = np.zeros((lx,ly))
time_map = np.asarray([i*tinc+tstart for i in range(nframes)])

# Useful functions for plotting cell fields
def iwrap(x):
    global lx
    remainder = x % lx
    if (remainder >= 0):
        return remainder
    return lx + remainder

def jwrap(y):
    global ly
    remainder = y % ly
    if (remainder >= 0):
        return remainder
    return ly + remainder

# Read position data
nlines = npoints + 2
reader = open(pos_file, 'r')

while True:
    # Read header section (including time info)
    for i in range(2):
        line = reader.readline()
    
    if (not line): break
    data = line.split()
    time = int(data[1])
    if (time > tend):
        break
    elif (time < tstart or (time-tstart)%tinc != 0):
        for i in range(npoints):
            line = reader.readline()
    else:
#        print("Reading data at timestep = {:d}".format(time))
        frame = (time-tstart)//tinc
        for n in range(npoints):
            line = reader.readline()
            data = line.split()
            
            # Read wrapped position data
            x = float(data[0])
            y = float(data[1])
            pos[frame][n] = (x,y)
            
reader.close()

mat = np.zeros((2,2))

# Read field data
def read_field(frame, field):
    global field_file_dir, field_file_name, tstart, tinc
    time = tstart+frame*tinc
    filename = field_file_dir + "/" + field_file_name + \
               (".dat.{:d}".format(time))
    with open(filename,'r') as field_reader:
        for line in field_reader:
            if (len(line) == 0 or line[0] == '#'): continue
            data = line.split()
            if (len(data) == 0): continue
            i = int(data[0])
            j = int(data[1])
            txx = float(data[xxcol])
            tyy = float(data[yycol])
            txy = float(data[xycol])
            # Make traceless tensor
            #htr = 0.5*(txx+tyy)
            #sxx = sign*(txx-htr)
            #sxy = sign*(txy-htr)
            # Compute principal eigenvalue and eigenvector
            #w = 0.5*(np.arctan2(-sxy,-sxx)+np.pi)
            #print(w)
            #px = np.cos(w)
            #py = np.sin(w)
            
            mat[0,0] = txx
            mat[0,1] = txy
            mat[1,0] = txy
            mat[1,1] = tyy
            
            w, v = np.linalg.eig(mat)
            k = 0 if w[0] > w[1] else 1
            vx[i,j] = 5.0*v[0,k] #2.0*px
            vy[i,j] = 5.0*v[1,k] #2.0*py
            field[i,j] = w[k] #np.sqrt(sxx*sxx+sxy*sxy)
        
# Make animation
# Plot settings
fontsize = 14 #20
#cmap = mplcm.viridis #mplcm.RdYlBu_r
#fcmap = mplcm.binary_r
cmap = mplcm.RdYlBu_r
fcmap = mplcm.RdYlBu_r
field_min = val_min
field_max = val_max
norm = mpl.colors.Normalize(vmin=val_min, vmax=val_max, clip=True)
mapper = mplcm.ScalarMappable(norm=norm, cmap=cmap)
mapper.set_array([])

fig, ax = plt.subplots()
ax.set_aspect('equal') # For maintaining aspect ratio
ax.set_xlim([0,lx])
ax.set_ylim([0,ly])

if (use_cbar):
    # 0.046, 0.04
#    cbar = plt.colorbar(mapper,fraction=0.15,pad=0.02,format="%.1e")
    cbar = plt.colorbar(mapper,fraction=0.15,pad=0.02) #,format="%.1e")
    cbar.set_ticks(np.arange(tic_start, tic_end+tic_inc/2.0, tic_inc))

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'FreeSans' # A font close to Helvetica

ax.tick_params(axis="both", labelsize=fontsize)
if (use_cbar):
    cbar.ax.tick_params(labelsize=fontsize)

# Draw borders but no axes and ticks
plt.tick_params(axis="both", which="both", bottom=False, top=False,
                labelbottom=False, right=False, left=False, labelleft=False)

# Set plot margins
#plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0, hspace=0)

# Get the artist for plotting the field
plt_field = ax.imshow(field.transpose(), origin="lower", vmin=field_min, 
                      vmax=field_max, cmap=fcmap, zorder=1)

# Get the artist for plotting centre of cells
plt_pts, = ax.plot(pos[frame,:,0], pos[frame,:,1], '.',
                   markersize=5, color="black", zorder=2)

# Get the artist for plotting polarity of cells
if (head):
    plt_vec = ax.quiver(xx, yy, vx, vy, units="xy", scale=5.0,
                        zorder=3, norm=norm, color="black")    
else:
    plt_vec = ax.quiver(xx, yy, vx, vy, units="xy", scale=5.0,
                        zorder=3, headlength=0, headaxislength=0, norm=norm, 
                        color="black")

# Get the artist for plotting the time label
if (use_label):
    plt_time_txt = ax.text(0.425,0.01, "", fontsize=14,
                           horizontalalignment="center",
                           transform=plt.gcf().transFigure)
    plt_time_txt.set_fontsize(fontsize)

def plot_data(frame):
    global pos, time_map
    
    # Plot the field
    read_field(frame,field)
    plt_field.set_data(field.transpose())

    # Plot centres of cells
    plt_pts.set_xdata(pos[frame,:,0])
    plt_pts.set_ydata(pos[frame,:,1])

    # Plot polarity of cells
    #plt_vec.set_offsets(vec[:,0:2])
    plt_vec.set_UVC(vx,vy)
    
    # Plot the time label
    if (use_label):
        plt_time_txt.set_text(r"$D_rt = {:.1f}$".format(time_map[frame]*Dr*dt))
    
    return plt_pts, plt_vec, plt_field

if (make_movie):
    if (print_to_screen):
        ani = FuncAnimation(fig, plot_data, np.arange(nframes), 
                            fargs=[], interval=1)
        plt.show()
    else:
        ani = FuncAnimation(fig, plot_data, np.arange(nframes), 
                            fargs=[], interval=1)
        Writer = animation.writers["ffmpeg"]
        #writer = Writer(fps=15, bitrate=1500)
        writer = Writer(fps=5, bitrate=1500)
        ani.save(out_file,writer=writer)
else:
    plot_data(0)
    if (print_to_screen):
        plt.show()
    else:
        plt.savefig(out_file, bbox_inches='tight', transparent=True)

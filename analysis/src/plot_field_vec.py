#!/bin/env python3
# plot_field_velocity.py

import sys

# This code must be run with python3!
if (sys.version_info < (3, 5)):
    print("This code must be run with Python version 3.5 or higher")
    sys.exit(1)

import numpy as np
#from skimage import measure
import matplotlib as mpl
import matplotlib.cm as mplcm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import rc
from matplotlib.animation import FuncAnimation
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

args = sys.argv
if (len(args) != 25):
    print("Usage: plot_field_vec.py npoints lx ly Dr dt vec_xcol vec_ycol vec_vcol vec_label vec_min vec_max tic_start tic_end tic_inc tstart tend tinc make_movie print_to_screen field_file_dir field_file_name pos_file polar_file out_file")
    sys.exit(1)

npoints = int(args.pop(1))
lx = int(args.pop(1))
ly = int(args.pop(1))
Dr = float(args.pop(1))
dt = float(args.pop(1))
vec_xcol = int(args.pop(1))
vec_ycol = int(args.pop(1))
vec_vcol = int(args.pop(1))
vec_label = args.pop(1)
vec_min = float(args.pop(1))
vec_max = float(args.pop(1))
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
polar_file = args.pop(1)
out_file = args.pop(1)
nframes = (tend-tstart)//tinc+1

use_label = 1 # 1
use_cbar = 1 # 1
head = 0 # 1
time_shift = 1 # 1 - report time values from the start time
vec_label = vec_label.replace("_"," ")

if (not make_movie):
    tend = tstart

# Data arrays
vec_len = (0.5*np.sqrt(lx*ly/npoints))
pos = np.zeros((nframes,npoints,2))
vec = np.zeros((nframes,npoints,5))
field = np.zeros((lx,ly))
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
        frame = (time-tstart)//tinc
        for n in range(npoints):
            line = reader.readline()
            data = line.split()
            
            # Read wrapped position data
            x = float(data[0])
            y = float(data[1])
            pos[frame][n] = (x,y)
            
reader.close()

# Read polar file
reader = open(polar_file, 'r')
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
            reader.readline()
    else:
        frame = (time-tstart)//tinc
        for n in range(npoints):
            line = reader.readline()
            data = line.split()
            xcm = pos[frame,n,0]
            ycm = pos[frame,n,1]
            vx = float(data[vec_xcol])
            vy = float(data[vec_ycol])
            v = np.sqrt(vx*vx+vy*vy)
            px = (vx/v)*vec_len
            py = (vy/v)*vec_len
            if (vec_vcol != -1):
                v = float(data[vec_vcol])
            vec[frame,n,0] = xcm-px/2.0
            vec[frame,n,1] = ycm-py/2.0
            vec[frame,n,2] = px
            vec[frame,n,3] = py
            vec[frame,n,4] = v
reader.close()
    
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
            val = float(data[2])
            field[i,j] = val
        
# Make animation
# Plot settings
fontsize = 14 #20
#cmap = mplcm.viridis #mplcm.RdYlBu_r
#fcmap = mplcm.binary_r
fcmap = mplcm.plasma
cmap = mplcm.Blues #mplcm.RdYlBu_r
field_min = 0.0 #2.6
field_max = 1.5 #3.2
norm = mpl.colors.Normalize(vmin=vec_min, vmax=vec_max, clip=True)
mapper = mplcm.ScalarMappable(norm=norm, cmap=cmap)
mapper.set_array([])

fig, ax = plt.subplots()
ax.set_aspect('equal') # For maintaining aspect ratio
ax.set_xlim([0,lx])
ax.set_ylim([0,ly])

if (use_cbar):
    cbar = plt.colorbar(mapper,fraction=0.15,pad=0.02) #,format="%.1e")
    cbar.set_ticks(np.arange(tic_start, tic_end+tic_inc/2.0, tic_inc))

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'FreeSans' # A font close to Helvetica

ax.tick_params(axis="both", labelsize=fontsize)
if (use_cbar):
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel(r"{:s}".format(vec_label), fontsize=fontsize, 
                       rotation=270)

# Draw borders but no axes and ticks
plt.tick_params(axis="both", which="both", bottom=False, top=False,
                labelbottom=False, right=False, left=False, labelleft=False)

# Set plot margins
#plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0, hspace=0)

# Get the artist for plotting the field
plt_field = ax.imshow(field.transpose(), origin="lower", vmin=field_min, 
                      vmax=field_max, cmap=fcmap, zorder=1)

# Get the artist for plotting polarity of cells
if (head):
    plt_vec = ax.quiver(vec[frame,:,0],vec[frame,:,1],vec[frame,:,2],
                        vec[frame,:,3],vec[frame,:,4],units="xy", zorder=3,
                        norm=norm, cmap=cmap)    
else:
    plt_vec = ax.quiver(vec[frame,:,0],vec[frame,:,1],vec[frame,:,2],
                        vec[frame,:,3],vec[frame,:,4],units="xy", zorder=3,
                        headlength=0, headaxislength=0, norm=norm, cmap=cmap)

# Get the artist for plotting the time label
if (use_label):
    plt_time_txt = ax.text(0.45,0.01, "", fontsize=fontsize,
                           horizontalalignment="center",
                           transform=plt.gcf().transFigure)

def plot_data(frame):
    global pos, time_map
    
    # Plot the field
    read_field(frame,field)
    plt_field.set_data(field.transpose())

    # Plot polarity of cells
    plt_vec.set_offsets(vec[frame,:,0:2])
    plt_vec.set_UVC(vec[frame,:,2],vec[frame,:,3],vec[frame,:,4])
    
    # Plot the time label
    if (use_label):
        if (time_shift):
            time = (time_map[frame]-tstart)*Dr*dt
        else:
            time = time_map[frame]*Dr*dt
        plt_time_txt.set_text(r"$D_rt = {:.1f}$".format(time))
    
    return plt_vec, plt_field

if (make_movie):
    if (print_to_screen):
        ani = FuncAnimation(fig, plot_data, np.arange(nframes), 
                            fargs=[], interval=1)
        plt.show()
    else:
        ani = FuncAnimation(fig, plot_data, np.arange(nframes), 
                            fargs=[], interval=1)
        Writer = animation.writers["ffmpeg"]
        writer = Writer(fps=10, bitrate=1000)
        ani.save(out_file,writer=writer)
else:
    plot_data(0)
    if (print_to_screen):
        plt.show()
    else:
        plt.savefig(out_file, bbox_inches='tight', transparent=True)

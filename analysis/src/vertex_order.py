# vertex_order.py
# A script to compute the p-fold complex shape order parameter for the cells
# in Voronoi representation

import sys

if (sys.version_info < (3, 5)):
    print("This code must be run with Python version 3.5 or higher")
    sys.exit(1)

import numpy as np
import scipy.spatial as spatial

args = sys.argv
if (len(args) != 12):
    print("Usage: vertex_order.py npoints lx ly xbuff ybuff degree",
          "tstart tend tinc pos_file out_file")
    sys.exit(1)

npoints = int(args.pop(1))
lx = int(args.pop(1))
ly = int(args.pop(1))
xbuff = float(args.pop(1))
ybuff = float(args.pop(1))
degree = int(args.pop(1))
tstart = int(args.pop(1))
tend = int(args.pop(1))
tinc = int(args.pop(1))
pos_file = args.pop(1)
out_file = args.pop(1)
xbuff *= lx
ybuff *= ly

periodic_loc = [(lx,-ly),(lx,0),(lx,ly),(0,-ly),(0,ly),(-lx,-ly),
                (-lx,0),(-lx,ly)]

lxmax = lx+xbuff
lxmin = -xbuff
lymax = ly+ybuff
lymin = -ybuff

pos = []

# Function to convert to complex exponential from Cartesian data for vertices
to_cmplx = np.vectorize(lambda x,y : complex(x,y)**degree)

# Read position data
with open(pos_file, 'r') as reader, open(out_file, 'w') as writer:
    while True:
        # Read header lines
        for i in range(2):
            line = reader.readline()
        if (not line): break
        data = line.split()
        time = int(data[1])

        if (time > tend): break
        if (time < tstart or (time-tstart) % tinc != 0):
            for i in range(npoints):
                reader.readline()
            continue

        # Clear position data
        pos.clear()
        
        # Add all points within the boundary first
        for n in range(npoints):
            line = reader.readline()
            data = line.split()
            x = float(data[0])
            y = float(data[1])
            pos.append((x,y))
        
        # Add points in the buffer region from the periodic image
        for n in range(npoints):
            for p in periodic_loc:
                xp = pos[n][0]+p[0]
                yp = pos[n][1]+p[1]
                if (xp < lxmax and x > lxmin and y < lymax and y > lymin):
                    pos.append((xp,yp))
        
        # Do Voronoi tesellation
        vor = spatial.Voronoi(pos)
        
        # Compute the complex order parameter
        writer.write("Cells: {:d}\n".format(npoints))
        writer.write("Timestep: {:d}\n".format(time))
        for n in range(npoints):
            region = vor.regions[vor.point_region[n]]
            if (-1 in region): continue # Ignore polygons extending to infinity

            # Compute CM
            pts = np.matrix([vor.vertices[i] for i in region])
            cm = pts.mean(0)

            # Shift points to CM frame
            pts -= cm

            # Compute normalisation factor
            norm = np.sum(np.linalg.norm(pts, ord=2, axis=1)**degree)
            order = np.sum(to_cmplx(pts[:,0],pts[:,1]))/norm
            writer.write("{:g} {:g}\n".format(np.abs(order), np.angle(order)))
            

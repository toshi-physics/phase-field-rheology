# init_poisson_disc.py
# Generate the coordinates of the cells' centres of mass randomly using the
# Poisson disc method, implementing the Bridson's algorithm

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("lx", type=int, help="size of lattice in x direction")
parser.add_argument("ly", type=int, help="size of lattice in y direction")
parser.add_argument("ncells", type=int, help="number of cells")
parser.add_argument("seed", type=int, help="seed for random generator")
parser.add_argument("outfile", help="name of output file")
args = parser.parse_args()

lx = args.lx
ly = args.ly
ncells = args.ncells
seed = args.seed
outfile = args.outfile
maxAttempt = 30
multiplier = 3.0

rng = np.random.RandomState(seed)

# Determine the minimum distance between cells
ln = int(np.ceil(np.sqrt(ncells*multiplier)))
d = min(lx/float(ln), ly/float(ln))
radius = 2**0.5*d
nx = int(lx/float(d))
ny = int(ly/float(d))
glx = d*nx
gly = d*ny
print(ln,d,radius,nx,ny,glx,gly)

# Position the cells at the centre of the box
xbuff = (lx-glx)*0.5
ybuff = (ly-gly)*0.5

# Create the grid of points
pos = [] # Actual positions of the cells
grid = [[None for j in range(ny)] for i in range(nx)]
active = [] # Active list of points for creating new points

# Helper functions
def dist(x1, y1, x2, y2):
    dx = x2-x1
    dy = y2-y1
    return np.sqrt(dx*dx+dy*dy)

def isValidPoint(x, y):
    if (x < 0.0 or x >= glx or y < 0.0 or y >= gly):
        return False
    
    # Check neighbouring cells (there is max one point per cell)
    xi = int(np.floor(x/d))
    yi = int(np.floor(y/d))
    i0 = max(xi-1, 0)
    i1 = min(xi+1, nx-1)
    j0 = max(yi-1, 0)
    j1 = min(yi+1, ny-1)
    for i in range(i0, i1+1):
        for j in range(j0, j1+1):
            if (grid[i][j] != None):
                px, py = grid[i][j]
                if (dist(x, y, px, py) < radius):
                    return False
    return True
    
def insertPoint(x, y):
    xi = int(np.floor(x/d))
    yi = int(np.floor(y/d))
    print(xi,yi,x,y)
    grid[xi][yi] = (x,y)


# Do sampling
x = rng.random()*glx
y = rng.random()*gly
insertPoint(x,y)
pos.append((x,y))
active.append((x,y))

while (len(active) > 0 and len(pos) < ncells):
    m = rng.randint(len(active))
    px, py = active[m]
    found = False
    for k in range(maxAttempt):
        theta = rng.random()*2.0*np.pi
        r = rng.random()*radius+radius
        newpx = px + r*np.cos(theta)
        newpy = py + r*np.sin(theta)
        if (not isValidPoint(newpx, newpy)):
            continue
        pos.append((newpx, newpy))
        insertPoint(newpx, newpy)
        active.append((newpx, newpy))
        found = True
        break
    if (not found):
        active.pop(m)

with open(outfile, "w") as writer:
    for i, xy in enumerate(pos):
        x, y = xy
        writer.write("{:d} {:.5f} {:.5f}\n".format(i, x+xbuff, y+ybuff))

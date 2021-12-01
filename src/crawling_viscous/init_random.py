# init_random.py
# Generate the coordinates of the cells' centres of mass randomly within
# a specified domain

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("lx", type=int, help="size of lattice in x direction")
parser.add_argument("ly", type=int, help="size of lattice in y direction")
parser.add_argument("ncells", type=int, help="number of cells")
parser.add_argument("radius", type=float, help="contact radius of each cell")
parser.add_argument("maxAttempt", type=int, help=
                    "max. number of attempts at generating a cell position")
parser.add_argument("seed", type=int, help="seed for random generator")
parser.add_argument("outfile", help="name of output file")
args = parser.parse_args()

lx = args.lx
ly = args.ly
ncells = args.ncells
rc = args.radius
maxAttempt = args.maxAttempt
seed = args.seed
outfile = args.outfile

rng = np.random.RandomState(seed)

pos = []

# Consider periodic boundaries
def sgnd(val):
    return (0.0 < val) - (val < 0.0)

def ddiff(width, d1, d2):
    dd1 = d1-d2
    add1 = np.abs(dd1)
    add2 = width-add1
    return dd1 if add1 < add2 else -sgnd(dd1)*add2

def rdiff(lx, ly, x1, x2, y1, y2):
    dx = ddiff(lx, x1, x2)
    dy = ddiff(ly, y1, y2)
    return np.sqrt(dx*dx+dy*dy)

reachedMaxAttempt = False
for i in range(ncells):
    print(i)
    if (reachedMaxAttempt):
        print("Reached max number of attempts - try using a different seed",
              "or a smaller radius")
        break
    
    nattempt = 0
    while (nattempt < maxAttempt):
        x = rng.random()*lx
        y = rng.random()*ly
        failed = False
        for px,py in pos:
            print(i, nattempt, x, y, px, py, rdiff(lx, ly, x, px, y, py))
            if (rdiff(lx, ly, x, px, y, py) < rc):
                failed = True
                break
        if (not failed):
            pos.append((x,y))
            break
        nattempt += 1
    if (nattempt == maxAttempt): 
        reachedMaxAttempt = True

with open (outfile, "w") as writer:
    for i,xy in enumerate(pos):
        writer.write("{:d} {:.5f} {:.5f}\n".format(i,*xy))

# plot_contour.py

import sys
import numpy as np
import cv2 as cv
from skimage import transform, measure
import matplotlib.pyplot as plt

args = sys.argv

if (len(args) != 4):
    print("Usage: plot_contour.py lx ly field_file")
    sys.exit(1)

lx = int(args.pop(1))
ly = int(args.pop(1))
field_file = args.pop(1)

phi = np.zeros((lx,ly))

def read_field(field, field_file):
    with open(field_file, 'r') as reader:
        for line in reader:
            data = line.split()
            if (len(data) == 0): continue
            i = int(data[0])
            j = int(data[1])
            field[i,j] = float(data[2])

read_field(phi, field_file)

maxphi = np.max(phi)
thres = 1.0
scale = 1.0
phi = transform.resize(phi, np.asarray(phi.shape)*scale)
contours = measure.find_contours(phi, 1.0)

fig, ax = plt.subplots()
#ax.imshow(phi)

for contour in contours:
    ax.plot(contour[:,1], contour[:,0])

plt.show()

#ret, phibin = cv.threshold(phi, thres, 1.0, cv.THRESH_BINARY)
#
#phibin = np.asarray(phibin, dtype=int)
#print(phibin)
#phibin *= 255
#contours,hierarchy = cv.findContours(phibin, cv.RETR_FLOODFILL, cv.CHAIN_APPROX_NONE)
#print(phi)
#print(phibin)
#print(contours)

#_,contours,_ = cv.findContours(phi, cv.RETR_TREE, )
#print(contours)

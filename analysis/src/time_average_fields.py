# time_average_fields.py

import sys
import math

args = sys.argv

if (len(args) != 8):
    print("Usage: time_average_fields.py ncols time_col tstart tend tinc",
          "data_file out_file")
    sys.exit(1)

ncols = int(args.pop(1))
time_col = int(args.pop(1))
tstart = int(args.pop(1))
tend = int(args.pop(1))
tinc = int(args.pop(1))
data_file = args.pop(1)
out_file = args.pop(1)

# Read data file and store average
avg = [0.0 for i in range(ncols)]
avgSq = [0.0 for i in range(ncols)]
n = 0

cols = [i for i in range(ncols)]
cols.pop(time_col)

with open(data_file, "r") as reader:
    for line in reader:
        # Ignore comments
        if (line.startswith("#")): continue
        data = line.split()
        
        # Ignore any empty lines
        if (data == []): continue
        
        time = int(data[time_col])
        
        # Ignore lines before tstart and after tend
        if (time > tend): break
        if (time < tstart or (time-tstart) % tinc != 0):  continue
        
        for i in cols:
            if (i != time_col):
                value = float(data[i])
                avg[i] += value
                avgSq[i] += value*value
        n += 1

n = float(n)
for i in cols:
    avg[i] /= max(n,1.0)
    avgSq[i] /= max(n,1.0)

# Compute error (use unbiased estimate for standard deviation)
with open(out_file, "w") as writer:
    if (n > 1.0):
        for i in cols:
            var = n / (n-1.0) * (avgSq[i] - avg[i]*avg[i])
            if (var < 0.0):
                print("Negative variance: var = {:.5f}".format(var))
            sigma = math.sqrt(abs(var)) 
            error = sigma / math.sqrt(n)
            writer.write("{:.5e} {:.5e} {:.5e} ".format(avg[i], sigma, error))
        writer.write("\n")
    else:
        for i in cols:
            writer.write("{:.5e} {:.5e} {:.5e}\n".format(avg[i], 0.0, 0.0))
        writer.write("\n")


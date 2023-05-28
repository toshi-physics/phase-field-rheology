# average_multi_files_fields.py
import sys
import math

args = sys.argv

if (len(args) < 6):
    print("Usage: average_multi_files.py ncols ref_col startpt",
          "output_file data_files")
    sys.exit(1)

args.pop(0) # Ignore self
ncols = int(args.pop(0))
ref_col = int(args.pop(0))
startpt = int(args.pop(0))
output_file = args.pop(0)

files = [open(i, "r") for i in args]
writer = open(output_file, "w")

cols = [i for i in range(ncols)]
has_ref_col = False
if (ref_col >= 0 and ref_col < ncols):
    cols.pop(ref_col)
    has_ref_col = True

for rows in zip(*files):
    ref = 0.0
    avg = [0.0 for i in range(ncols)]
    avgSq = [0.0 for i in range(ncols)]
    dataOK = False
    n = 0

    for j,line in enumerate(rows):
        if (line.startswith("#")): continue
        data = line.split()
        if (data == []):  continue
        # Convert all columns to numerical values
        for i in range(ncols):
            data[i] = float(data[i])
            
        if (has_ref_col):
            ref = data[ref_col]

        if ((not has_ref_col) or ref >= startpt):
            for i in cols:
                value = data[i]
                avg[i] += value
                avgSq[i] += value*value
            n += 1
            dataOK = True

    if (dataOK):
        n = float(n)
        if (n > 1.0):
            if (has_ref_col):
                writer.write("{:f} ".format(ref))
            for i in cols:
                avg[i] /= n
                avgSq[i] /= n
                
                # Use un-biased estimate of variance
                var = n / (n-1.0) * (avgSq[i] - avg[i]*avg[i])
                if (var < 0.0):
                    print("Negative variance: var = {:e}".format(var))
                sigma = math.sqrt(abs(var))
                error = sigma / math.sqrt(n)
                output = "{:e} {:e} {:e} ".format(avg[i], sigma, error)
                writer.write(output)
            writer.write("\n")
        else:
            for i in cols:
                output = "{:e} {:e} {:e} ".format(avg[i], 0.0, 0.0)
                writer.write(output)
            writer.write("\n")

for f in files:
    f.close()
writer.close()

# average_multi_files.py
import sys
import math

args = sys.argv

if (len(args) < 7):
    print("Usage: average_multi_files.py ref_col avg_col err_col",
          "startpt output_file data_files")
    sys.exit(1)

args.pop(0) # Ignore self
ref_col = int(args.pop(0))
avg_col = int(args.pop(0))
err_col = int(args.pop(0))
startpt = int(args.pop(0))
output_file = args.pop(0)

files = [open(i, "r") for i in args]
writer = open(output_file, "w")

has_ref_col = False
if (ref_col >= 0):
    has_ref_col = True
has_err_col = False
if (err_col >= 0):
    has_err_col = True

for rows in zip(*files):
    ref = 0.0
    avg = 0.0
    avgSq = 0.0
    error = 0.0
    n = 0
    hasData = False

    for line in rows:
        if (line.startswith("#")): continue
        data = line.split()
        if (data == []): continue # Ignore any lines start with \n
            
        if (has_ref_col):
            ref = float(data[ref_col])
        
        if ((not has_ref_col) or ref >= startpt):
            value = float(data[avg_col])
            avg += value
            avgSq += value*value
            n += 1
        
        if (has_err_col):
            sigma = float(data[err_col])
            error += sigma*sigma

        hasData = True

    if (hasData):
        n = float(n)
        avg /= n
        avgSq /= n

        # use un-biased estimate of variance
        if (n > 1.0):
            var = n / (n-1.0) * (avgSq - avg*avg)
            if (var < 0.0):
                print("Negative variance: var = {:e}".format(var))
            sigma = math.sqrt(abs(var)) 
        else:
            var = 0
            sigma = 0
        
        if (has_err_col):
            error = math.sqrt(error) / n
        else:
            error = sigma / math.sqrt(n)

        if (has_ref_col):
            output = "{:e} {:e} {:e} {:e}\n".format(ref, avg, sigma, error)
        else:
            output = "{:e} {:e} {:e}\n".format(avg, sigma,  error)
        writer.write(output)

for f in files:
    f.close()
writer.close()

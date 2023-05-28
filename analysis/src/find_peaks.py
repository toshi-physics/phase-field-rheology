# find_peaks.py
# A script to find local maxima in a curve

import sys

args = sys.argv

if (len(args) != 5):
    print("Usage: find_peaks.py ref_col val_col data_file out_file")
    sys.exit(1)

ref_col = int(args.pop(1))
val_col = int(args.pop(1))
data_file = args.pop(1)
out_file = args.pop(1)

diff = 0
has_prev_val = False
has_prev_diff = False
prev_val = 0.0
prev_diff = 0.0
prev_ref = 0.0

with open(data_file, 'r') as reader, open(out_file, 'w') as writer:
    for line in reader:
        if (line.startswith("#")): continue
        data = line.split()
        if (len(data) == 0): continue
        
        ref = float(data[ref_col])
        val = float(data[val_col])

        if (has_prev_val):
            diff = (val - prev_val)
            if (not has_prev_diff):
                has_prev_diff = True
            elif (prev_diff > 0.0 and diff < 0.0): 
                writer.write("{:f} {:f}\n".format(prev_ref, prev_val))
            prev_diff = diff
        else:
            has_prev_val = True
        prev_val = val
        prev_ref = ref
        

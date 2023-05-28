# spearman_correlation.py
# A program which computes the Spearman correlation coefficient of two datasets

import sys
import scipy
import scipy.stats
import math

args = sys.argv

if (len(args) != 9):
    print("Usage: spearman_correlation.py size index_col1 value_col1",
          "index_col2 value_col2 array_file1 array_file2 out_file")
    sys.exit(1)

size = int(args.pop(1))
index_col1 = int(args.pop(1))
value_col1 = int(args.pop(1))
index_col2 = int(args.pop(1))
value_col2 = int(args.pop(1))
array_file1 = args.pop(1)
array_file2 = args.pop(1)
out_file = args.pop(1)

def read_array(filename, array, index_col, value_col):
    index = 0
    with open(filename, 'r') as f:
        for line in f:
            if (line.startswith("#")): continue
            data = line.split()
            if (data == []): continue
            if (index_col != -1):
                index = int(data[index_col])
            value = float(data[value_col])
            array[index] = value
            if (index_col == -1):
                index += 1

array1 = [0.0 for i in range(size)]
array2 = [0.0 for i in range(size)]

# Read in arrays 
read_array(array_file1, array1, index_col1, value_col1)
read_array(array_file2, array2, index_col2, value_col2)

# Remove nan or inf entries
arrays = [[a,b] for a,b in zip(array1,array2) if not 
          (math.isnan(a) or math.isinf(a) or math.isnan(b) or math.isinf(b))]
array1 = [x[0] for x in arrays]
array2 = [x[1] for x in arrays]

array1 = scipy.array(array1)
array2 = scipy.array(array2)

# Compute correlation
corr = scipy.stats.spearmanr(array1, array2)

with open(out_file, 'w') as writer:
    writer.write("{:.5g} {:.5g}\n".format(*corr))



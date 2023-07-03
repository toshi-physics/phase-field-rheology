#!/bin/bash

# Directory of this script
sh_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

if (( $# != 10 )); then
    echo "Usage: phase_field_multi_init.sh d_start d_end d_inc pe_start pe_end pe_inc run_start run_end run_inc dir"
    exit 1
fi

d_start=$1
d_end=$2
d_inc=$3
pe_start=$4
pe_end=$5
pe_inc=$6
run_start=$7
run_end=$8
run_inc=$9
dir=${10}

pyx=python3

d=$($pyx -c "print('{:.3f}'.format($d_start))")
pe=$($pyx -c "print('{:.3f}'.format($pe_start))")
run=$run_start

while (( $(bc <<< "$d <= $d_end") ))
do
    pe=$($pyx -c "print('{:.3f}'.format($pe_start))")
    while (( $(bc <<< "$pe <= $pe_end") ))
    do
	for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	do
	    echo "Creating files for d = $d pe = $pe run = $run"
	    ./phase_field_init.sh $d $pe $a $run $dir
	done
	pe=$($pyx -c "print('{:.3f}'.format($pe+$pe_inc))")
    done
    d=$($pyx -c "print('{:.3f}'.format($d+$d_inc))")
done

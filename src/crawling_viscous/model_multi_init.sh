#!/bin/bash

# Directory of this script
sh_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

if (( $# != 13 )); then
    echo "Usage: model_multi_init.sh d_start d_end d_inc pe_start pe_end pe_inc eta_start eta_end eta_inc run_start run_end run_inc dir"
    exit 1
fi

d_start=$1
d_end=$2
d_inc=$3
pe_start=$4
pe_end=$5
pe_inc=$6
eta_start=$7
eta_end=$8
eta_inc=$9
run_start=${10}
run_end=${11}
run_inc=${12}
dir=${13}

pyx=python3

d=$($pyx -c "print('{:.2f}'.format($d_start))")
pe=$($pyx -c "print('{:.2f}'.format($pe_start))")
eta=$($pyx -c "print('{:.2f}'.format($eta_start))")
run=$run_start

while (( $(bc <<< "$d <= $d_end") ))
do
    pe=$($pyx -c "print('{:.2f}'.format($pe_start))")
    while (( $(bc <<< "$pe <= $pe_end") ))
    do
	eta=$($pyx -c "print('{:.2f}'.format($eta_start))")
	while (( $(bc <<< "$eta <= $eta_end") ))
	do
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		echo "Creating files for d = $d Pe = $pe eta = $eta run = $run"
		./model_init.sh $d $pe $eta $run $dir
	    done
	    eta=$($pyx -c "print('{:.2f}'.format($eta+$eta_inc))")
	done
	pe=$($pyx -c "print('{:.2f}'.format($pe+$pe_inc))")
    done
    d=$($pyx -c "print('{:.2f}'.format($d+$d_inc))")
done

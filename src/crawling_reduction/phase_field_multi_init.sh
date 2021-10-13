#!/bin/bash

if [[ $# != 13 ]]; then
    echo "Usage: phase_field_multi_init.sh d_start d_end d_inc pe_start pe_end pe_inc a_start a_end a_inc run_start run_end run_inc dir"
    exit 1
fi

d_start=$1
d_end=$2
d_inc=$3
pe_start=$4
pe_end=$5
pe_inc=$6
a_start=$7
a_end=$8
a_inc=$9
run_start=${10}
run_end=${11}
run_inc=${12}
dir=${13}

pyx=python3

d=$($pyx -c "print('{:.3f}'.format($d_start))")
pe=$($pyx -c "print('{:.3f}'.format($pe_start))")
a=$($pyx -c "print('{:.3f}'.format($a_start))")
run=$run_start

while (( $(bc <<< "$d <= $d_end") ))
do
    pe=$($pyx -c "print('{:.3f}'.format($pe_start))")
    while (( $(bc <<< "$pe <= $pe_end") ))
    do
	a=$($pyx -c "print('{:.3f}'.format($a_start))")
	while (( $(bc <<< "$a <= $a_end") ))
	do
	    for (( run=$run_start; $run<=$run_end; run+=$run_inc ))
	    do
		echo "Creating files for d = $d pe = $pe a = $a run = $run"
		./phase_field_init.sh $d $pe $a $run $dir
	    done
	    a=$($pyx -c "print('{:.3f}'.format($a+$a_inc))")
	done
	pe=$($pyx -c "print('{:.3f}'.format($pe+$pe_inc))")
    done
    d=$($pyx -c "print('{:.3f}'.format($d+$d_inc))")
done

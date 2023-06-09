#!/bin/bash
# phase_field_init.sh

# Directory of this script
sh_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

if (( $# != 4 )); then
    echo "Usage: phase_field_init.sh deform peclet run run_dir"
    exit 1
fi

deform=$1    # Deformability (d = epsilon/alpha)
peclet=$2    # Peclet number (Pe = v/(R*Dr))
run=$3       # Trial number
run_dir=$4   # Run directory

pyx=python3 # Choose between python2 and python3

# Format input params
deform=$($pyx -c "print('{:.3f}'.format($deform))")
peclet=$($pyx -c "print('{:.3f}'.format($peclet))")

# A function for generating random numbers
max_seed=1000000
function get_rand(){
    # Generate a 4-byte random integer using urandom
    rand=$(od -vAn -N4 -tu4 < /dev/urandom)
    echo $rand
}

# Set the model parameters
ncells=100 # 100
ncell_x=10 # 10
ncell_y=10 # 10
confine_radius=8.0 # 8.0
init_radius=7.0 # 7.0
ideal_radius=12.0 # 12.0
cell_lx=41 # 41
cell_ly=41 # 41
phi0=1.0 # 1.0
lambda=1.0 # 1.0 (rougly 1800 in the old model)
epsilon=0.1 # 0.1
omega=0.1 # 0.1
kappa_active=0.1 # 0.1 active shear stress
rotate_diff=0.0001 # 0.0001
friction=1.0 # 1.0
mobility=0.1 # 0.1
overlap=1 # Do overlap calculation or not

nsteps=21000000 # 20000000 # Add another 10^6 steps for equilibration
nequil=10000 # 10000
delta_t=0.5 # 0.5
dump_cm_freq=1000 # 1000
dump_bulk_cm_freq=1000 # 1000
dump_gyr_freq=1000 # 1000
dump_gyr_field_freq=10000 # 1000
dump_vel_freq=1000 # 1000
dump_vel_field_freq=10000 # 1000
dump_deform_freq=1000 # 1000
dump_deform_field_freq=10000 # 1000
dump_field_freq=10000 # 100000
dump_cell_field_freq=10000 # 100000
dump_index_field_freq=10000 # 100000
dump_shape_freq=1000 # 1000
dump_neighbour_freq=1000 # 1000
dump_energy_freq=1000 # 1000
dump_overlap_freq=1000 # 1000
equildump_cm_freq=1000 # 1000
equildump_gyr_freq=10000 # 10000
equildump_field_freq=10000 # 10000
seed=$(get_rand)

# Set alpha based on deformability
alpha=$($pyx -c "print('{:f}'.format(1.5*$epsilon/$deform))")
kappa=$($pyx -c "print('{:f}'.format($alpha*2.0))")

# Set motility based on Peclet number, rotatioal diff, and ideal radius
motility=$($pyx -c "print('{:f}'.format($peclet*$rotate_diff*$ideal_radius))")

# Set a hexagonal lattice
tmp_cm_file="cm_$seed.tmp"
tmp_shape_file="shape_$seed.tmp"
#size=$($pyx triangle.py $ncell_x $ncell_y $confine_radius $tmp_cm_file)
#size=$($pyx square.py $ncell_x $ncell_y $confine_radius $tmp_cm_file)
#size=$($pyx square.2.py 160 138 $ncell_x $ncell_y $confine_radius $tmp_cm_file)
cutoff=0.2
seed_2=$(get_rand)
size=$($pyx triangle_noise.py $ncell_x $ncell_y $confine_radius $cutoff $seed_2 $tmp_cm_file)
#size=$(python triangle_stretched_noise.py $ncell_x $ncell_y $confine_radius $cutoff $seed_2 $tmp_cm_file)

lx=$(echo $size | awk '{print $3}')
ly=$(echo $size | awk '{print $6}')
python cell_shape.py $cell_lx $cell_ly $phi0 $tmp_shape_file circle $init_radius

# Create run directory and set file names
sim_name="cell_N_${ncells}_d_${deform}_Pe_${peclet}_run_${run}"
run_dir="${run_dir}/${sim_name}/"

if [ ! -d $run_dir ]; then
    mkdir -p $run_dir
fi

cm_file="cm_${sim_name}.in"
shape_file="shape_${sim_name}.in"
params_file="params_${sim_name}.txt"
#equildump_cm_file="pos-equil_${sim_name}.dat"
#equildump_gyr_file="gyr-equil_${sim_name}.dat"
#equildump_field_file="field-equil_${sim_name}.dat"
dump_cm_file="pos_${sim_name}.dat"
dump_gyr_file="gyr_${sim_name}.dat"
dump_gyr_field_file="gyr-field_${sim_name}.dat"
dump_vel_file="vel_${sim_name}.dat"
dump_vel_field_file="vel-field_${sim_name}.dat"
dump_deform_file="deform_${sim_name}.dat"
dump_deform_field_file="deform-field_${sim_name}.dat"
dump_field_file="field_${sim_name}.dat"
dump_cell_field_file="cell-field_${sim_name}"
dump_index_field_file="index-field_${sim_name}.dat"
dump_bulk_cm_file="pos-bulk_${sim_name}.dat"
dump_shape_file="shape_${sim_name}.dat"
dump_neighbour_file="neigh_${sim_name}.dat"
dump_energy_file="energy_${sim_name}.dat"
dump_overlap_file="olap_${sim_name}.dat"

# Copy the template file
params_file=${run_dir}/$params_file
mv $tmp_cm_file ${run_dir}/$cm_file
mv $tmp_shape_file ${run_dir}/$shape_file

# Replace macros in template with input values
echo \
"phi0 = ${phi0}
cellRadius = ${ideal_radius}
cahnHilliardCoeff = ${alpha}
surfaceTensionCoeff = ${kappa}
volumePenaltyCoeff = ${lambda}
repulsionCoeff = ${epsilon}
adhesionCoeff = ${omega}
frictionCoeff = ${friction}
activeShearCoeff = ${kappa_active}
diffusionCoeff = ${rotate_diff}
cellLx = ${cell_lx}
cellLy = ${cell_ly}
lx = ${lx}
ly = ${ly}
nsteps = ${nsteps}
nequil = ${nequil}
ncells = ${ncells}
dt = ${delta_t}
mobility = ${mobility}
motility = ${motility}
cm_file = ${cm_file}
shape_file = ${shape_file}
seed = ${seed}
overlap = ${overlap}
" > $params_file

# Set dumps
function add_dump() {
    params=$1; file=$2;
    if [ "$file" ]; then
	echo "$params $file" >> $params_file 
    fi
}
#add_dump "dump_cm $equildump_cm_freq 0 equil" $equildump_cm_file
#add_dump "dump_gyr $equildump_gyr_freq 0 equil" $equildump_gyr_file
#add_dump "dump_field $equildump_field_freq 0 equil" $equildump_field_file
add_dump "dump_cm $dump_cm_freq 0 main" $dump_cm_file
#add_dump "dump_gyr $dump_gyr_freq 0 main" $dump_gyr_file
#add_dump "dump_gyr_field $dump_gyr_field_freq 0 main" $dump_gyr_field_file
#add_dump "dump_vel $dump_vel_freq 0 main" $dump_vel_file
#add_dump "dump_vel_field $dump_vel_field_freq 0 main" $dump_vel_field_file
#add_dump "dump_deform $dump_deform_freq 0 main" $dump_deform_file
#add_dump "dump_deform_field $dump_deform_field_freq 0 main" $dump_deform_field_file
add_dump "dump_field $dump_field_freq 0 main" $dump_field_file
#add_dump "dump_index_field $dump_index_field_freq 0 main" $dump_index_field_file
add_dump "dump_bulk_cm $dump_bulk_cm_freq main" $dump_bulk_cm_file
add_dump "dump_shape 4 25 4.0 5 31 $dump_shape_freq 0 main" $dump_shape_file
#add_dump "dump_neighbour $dump_neighbour_freq 0 main" $dump_neighbour_file
#add_dump "dump_energy $dump_energy_freq 0 main" $dump_energy_file
#add_dump "dump_overlap $dump_overlap_freq 0 main" $dump_overlap_file

#for (( i=0; $i<$ncells; i++ ))
#do
#    add_dump "dump_cell_field $i $dump_cell_field_freq 0 main" "${dump_cell_field_file}_cell_${i}.dat"
#done

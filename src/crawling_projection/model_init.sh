#!/bin/bash
# phase_field_init.sh

# Directory of this script
sh_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"

if (( $# != 5 )); then
    echo "Usage: model_init.sh deform peclet eta0 run run_dir"
    exit 1
fi

deform=$1    # Deformability (d = epsilon/alpha)
peclet=$2    # Peclet number (Pe = v/(R*Dr))
eta0=$3       # \eta_0 (\eta_0 = \eta / (\gamma_v * R^2))
run=$4       # Trial number
run_dir=$5   # Run directory

pyx=python3 # Choose between python2 and python3

# Required scripts and programs
init_rand_py="${sh_dir}/init_random.py"

# Format input params
deform=$($pyx -c "print('{:.2f}'.format($deform))")
peclet=$($pyx -c "print('{:.2f}'.format($peclet))")
eta0=$($pyx -c "print('{:.2f}'.format($eta0))")

# A function for generating random numbers
max_seed=1000000
function get_rand(){
    # Generate a 4-byte random integer using urandom
    rand=$(od -vAn -N4 -tu4 < /dev/urandom)
    echo $rand
}

# Set the model parameters
ncells=100 # 300
init_radius=5.0
ideal_radius=8.0 # 12.0
lx=145
ly=145 # 435
cell_lx=35 # 41
cell_ly=35 # 41
max_cell_lx=35 # 41
max_cell_ly=35 # 41
phi0=1.0 # 1.0
lambda=0.1 # 1.0 (rougly 1800 in the old model)
epsilon=0.01 # 0.1
omega=0.0 # 0.1
xi=1.0 # 1.0 thickness
kappa_active=0.0 # 0.1 active shear stress
#eta=1.0 # first coefficient of viscosity
#zeta=1.0 # second coefficient of viscosity
rotate_diff=0.001 # 0.0001
friction=$($pyx -c "print('{:f}'.format(1.0/$ideal_radius**2))")
mobility=1.0 # 1.0 0.1
eta=$($pyx -c "print('{:f}'.format($eta0/($ideal_radius**2*$friction)))")
zeta=$eta
overlap=1 # Do overlap calculation or not

nsteps=200000 #11000000 # 20000000 # Add another 10^6 steps for equilibration
nshear=20001 #same as equilibration time rn 
nequil=10000 # 100000
delta_t=1.0 # 0.5
dump_cm_freq=100 #1000 # 1000
dump_bulk_cm_freq=100 #1000 # 1000
dump_gyr_freq=100 #1000 # 1000
dump_gyr_field_freq=1000 #10000 # 1000
dump_vel_freq=100 #1000 # 1000
dump_vel_field_freq=1000 #10000 # 1000
dump_deform_freq=100 #1000 # 1000
dump_deform_field_freq=1000 #10000 # 1000
dump_field_freq=1000 #10000 # 100000
dump_cell_field_freq=1000 #10000 # 100000
dump_index_field_freq=1000 #10000 # 100000
dump_shape_freq=100 #1000 # 1000
dump_neighbour_freq=100 #1000 # 1000
dump_energy_freq=100 #1000 # 1000
dump_overlap_freq=100 #1000 # 1000
equildump_cm_freq=100 #1000 # 1000
equildump_gyr_freq=10000 # 10000
equildump_field_freq=10000 # 10000
seed=$(get_rand)

# Set kappa based on deformability
kappa=$($pyx -c "print('{:f}'.format(3.0*$epsilon/$deform))")

# Set motility based on Peclet number, rotatioal diff, and ideal radius
motility=$($pyx -c "print('{:f}'.format($peclet*$rotate_diff*$ideal_radius))")

# Set shearrate based on motility, comparable for now?
shearrate=$($pyx -c "print('{:f}'.format($motility))")

# Create run directory and set file names
sim_name="cell_N_${ncells}_d_${deform}_Pe_${peclet}_e_${eta0}_run_${run}"
run_dir="${run_dir}/${sim_name}/"

if [ ! -d $run_dir ]; then
    mkdir -p $run_dir
fi

cell_file="${sim_name}.in"
cell_file_path="${run_dir}/${cell_file}"

# Generate random positions for cells
seed_2=$(get_rand)
max_attempt=10000
$pyx $init_rand_py $lx $ly $ncells $init_radius $max_attempt $seed_2 $cell_file_path

params_file="${run_dir}/params_${sim_name}.txt"
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

# Replace macros in template with input values
echo \
"cellRadius = ${ideal_radius}
thickness = ${xi}
cahnHilliardCoeff = ${kappa}
volumePenaltyCoeff = ${lambda}
repulsionCoeff = ${epsilon}
adhesionCoeff = ${omega}
frictionCoeff = ${friction}
activeShearCoeff = ${kappa_active}
firstViscousCoeff = ${eta}
secondViscousCoeff = ${zeta}
diffusionCoeff = ${rotate_diff}
cellLx = ${cell_lx}
cellLy = ${cell_ly}
lx = ${lx}
ly = ${ly}
nsteps = ${nsteps}
nshear = ${nshear}
nequil = ${nequil}
ncells = ${ncells}
dt = ${delta_t}
mobility = ${mobility}
motility = ${motility}
shearrate = ${shearrate}
cell_file = ${cell_file}
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
add_dump "dump_bulk_cm $dump_bulk_cm_freq main" $dump_bulk_cm_file
add_dump "dump_gyr $dump_gyr_freq 0 main" $dump_gyr_file
#add_dump "dump_gyr_field $dump_gyr_field_freq 0 main" $dump_gyr_field_file
add_dump "dump_vel $dump_vel_freq 0 main" $dump_vel_file
#add_dump "dump_vel_field $dump_vel_field_freq 0 main" $dump_vel_field_file
add_dump "dump_deform $dump_deform_freq 0 main" $dump_deform_file
add_dump "dump_deform_field $dump_deform_field_freq 0 main" $dump_deform_field_file
add_dump "dump_field $dump_field_freq 0 main" $dump_field_file
#add_dump "dump_index_field $dump_index_field_freq 0 main" $dump_index_field_file
add_dump "dump_shape 4 25 4.0 5 31 $max_cell_lx $max_cell_ly $dump_shape_freq 0 main" $dump_shape_file
#add_dump "dump_neighbour $dump_neighbour_freq 0 main" $dump_neighbour_file
#add_dump "dump_energy $dump_energy_freq main" $dump_energy_file
#add_dump "dump_overlap $dump_overlap_freq 0 main" $dump_overlap_file

#for (( i=0; $i<$ncells; i++ ))
#do
#    add_dump "dump_cell_field $i $dump_cell_field_freq 0 main" "${dump_cell_field_file}_cell_${i}.dat"
#done

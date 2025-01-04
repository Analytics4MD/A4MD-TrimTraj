#!/bin/bash

# Parameters
DEBUG="False"  # Debug mode
current_datetime=$(date +"%Y-%m-%d_%H-%M-%S")

# Check if the results folder exists
if [ ! -d "./results" ]; then
    # Create the results folder
    mkdir "./results"
fi

termination_criterion="termination_ess"  # "termination_lev" or "termination_ess"
stride=1  # Stride for loading trajectories
sampling_rate=0.002  # In ps

# Trajectories files
input_dir=("./files/xtc_files/output.checkpoint" "./files/xtc_files/output")
input_dirs=$(IFS=,; echo "${input_dir[*]}")

top_file="./files/xtc_files/boxed.pdb"  # Topology file for the trajectories


# LEVs file
levs_df="./files/xtc_files/levs.trajs40_stride1.pkl"

# LEV
range_min=4000  # LEV range
range_max=4400  # LEV range
range_tolerance=5  # LEV tolerance to remain in range (in %)

# declare -a window_lev
window_lev=(200) # Sliding window size for LEV
stable_th_lev=(80)  # Minimum number of frames considered in window to decide termination for LEV
# if multiple terminations are required use the following format
# window_lev=(200 250 300) # Sliding window size for LEV
# stable_th_lev=(80 90 100)  # Minimum number of frames considered in window to decide termination for LEV
window_levs=$(IFS=,; echo "${window_lev[*]}")
stable_th_levs=$(IFS=,; echo "${stable_th_lev[*]}")

# ESS
## TODO: add example for multiple values
var_th_ess=(10) # Minimum average differential variability threshold for ESS (in %)
window_ess=(250) # Sliding window size for ESS
window_esss=$(IFS=,; echo "${window_ess[*]}")
var_th_esss=$(IFS=,; echo "${var_th_ess[*]}")


# #  ----- Download Full trajectories -----
# Execute download_data.py
if [ ! -d "./files/xtc_files" ]; then
    # Create the annotations folder
    echo "Downloading data..."
    python download_data.py
    echo "Data downloaded successfully!"
fi
#  ----- Generate Annotations -----

# Check if the annotations folder exists
if [ ! -d "./results/annotations" ]; then
    # Create the annotations folder
    mkdir "./results/annotations"
fi
script_name="annotations_${current_datetime}"
annotations_output_folder="./results/annotations/${script_name}"
mkdir -p "$annotations_output_folder"

# Execute generate_annotations_def.py
# python generate_annotations_def.py --termination_criterion "$termination_criterion" --DEBUG "$DEBUG"
echo "Generating annotations..."
# python generate_annotations_def.py $termination_criterion $DEBUG $stride $sampling_rate $range_min $range_max $range_tolerance $annotations_output_folder $input_dirs_1 $input_dirs_2 $levs_df $window_levs $stable_th_levs $var_th_esss $window_esss 
python generate_annotations_def.py $termination_criterion $DEBUG $stride $sampling_rate $range_min $range_max $range_tolerance $annotations_output_folder $input_dirs $levs_df $window_levs $stable_th_levs $var_th_esss $window_esss 
echo "Annotations generated successfully!"
echo "Output folder: $annotations_output_folder"

#  ----- Generate Trajectories -----
# # Definitions
first_batch=0  # FS1 or FS2 runs first (ie 1 or 2 are trimmed), 0 means trim everything, -1 means dont trim
trim_first_rand="False"  # Terminate at a random point
trial=0  # Identifier of random trial [0, ..., 4]
run_last_full="False"  # Let restarted trajectory run to completion
    

# Check if the trajectories folder exists
if [ ! -d "./results/trajectories" ]; then
    # Create the trajectories folder
    mkdir "./results/trajectories"
fi
script_name="trajectories_${current_datetime}"
trajectories_output_folder="./results/trajectories/${script_name}"
mkdir -p "$trajectories_output_folder"

# Execute create_trajectory_files_def.py
echo "Generating trajectories..."
# python create_trajectory_files_def.py $termination_criterion $DEBUG $trajectories_output_folder $stride $top_file $input_dirs_1 $input_dirs_2 $annotations_output_folder $first_batch $trim_first_rand $trial $run_last_full
python create_trajectory_files_def.py $termination_criterion $DEBUG $trajectories_output_folder $stride $top_file $input_dirs $annotations_output_folder $first_batch $trim_first_rand $trial $run_last_full
echo "Trajectories generated successfully!"


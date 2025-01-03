from datetime import datetime
import math
import os
import sys

import mdtraj as md
import numpy as np
import pandas as pd
import multiprocessing


from loguru import logger



## -------------------- Functions -------------------- ##


def load_trajectory(input_dirs, filename):
    trajectories = []

    # Prepare trajectory path
    traj_path = os.path.join(input_dirs, filename)
    tokens = filename.split("_")

    if len(tokens) != 4:
        return

    conformation = tokens[1]
    traj_id = tokens[3]

    if int(traj_id) > 20:  # Skipping unstable trajectories for now
        return

    if not os.path.exists(traj_path + "/out_md/traj_comp_whole.xtc"):
        return

    # n_trajectories += 1

    # Load trajectory
    t = md.load(
        traj_path + "/out_md/traj_comp_whole.xtc",
        stride=stride,
        top=traj_path + "/setup/conf.pdb",
    )
    logger.debug(
        "Loaded trajectory {} in subset {} from conformation {} with stride {}: {} frames ({} ns)".format(
            traj_id,
            input_dirs,
            conformation,
            stride,
            len(t),
            len(t) * to_ns,
        )
    )
    trajectories.append(
        [
            traj_id,
            input_dirs,
            conformation,
            stride,
            len(t),
            len(t) * to_ns,
            t,
        ]
    )

    columns = [
        "traj_id",
        "subset",
        "conformation",
        "stride",
        "n_frames",
        "length",
        "mdtraj",
    ]
    trajectory_df = pd.DataFrame(trajectories, columns=columns)

    # logger.debug("{} trajectories".format(len(trajectories)))
    return trajectory_df


def n_effective_FTZ(x):
    n = len(x)
    xbar = np.mean(x)
    xnorm = x - xbar
    r = []
    nc = 0
    neff = n

    for k in range(n - 1):
        r.append(sum(xnorm[0 : n - 1 - k] * xnorm[k + 1 : n]) / sum(xnorm * xnorm))
        if r[k] < 0:
            nc = k - 1
            break

    ss = 0
    for k in range(nc):
        ss = ss + r[k]

    neff = n / (1 + 2 * ss)
    neff = min(neff, n)
    return neff


def analyse_trajectory(traj, evs,stable_th_lev,window_lev,tt, window_ess, var_th_ess):
    end_ev = traj.n_frames
    stable_evs = [0] * window_lev  # Circular buffer for marking stable frames
    last = 0

    end_ess = traj.n_frames
    ess_by_frame = []
    precision_by_frame = []
    x = []

    end_both = traj.n_frames
    ev = False
    ess = False

    # Study latest eigenvalues
    # for frame in range(5):
    for frame in range(traj.n_frames):
        # Calculate ESS
        ess_by_frame.append(
            n_effective_FTZ(evs[max(0, frame + 1 - window_ess) : frame + 1])
        )
        x.append(frame)

        # Calculate precision of ESS
        # precision_by_frame.append( math.sqrt( ess_by_frame[frame] )  / np.std(evs[max(0,frame+1-window):frame+1]))

        # Check termination conditions based on variation of ESS
        if end_ess == traj.n_frames and frame + 1 >= window_ess and ess is False:
            dx = np.diff(x[-window_ess:])
            dy = np.diff(ess_by_frame[-window_ess:])

            d = dy / dx
            davg = np.average(d)
            if davg > var_th_ess / 100:
                end_ess = frame + 1
                ess = True

        # Check termination conditions based on eigenvalue range
        if end_ev == traj.n_frames and frame + 1 >= window_lev and ev is False:
            if (
                evs[frame] >= range_min - range_th
                and evs[frame] <= range_max + range_th
            ):
                stable_evs[last] = 1
            else:
                stable_evs[last] = 0

            if np.sum(stable_evs) >= (stable_th_lev / 100) * window_lev:
                if tt == 0 or tt == 20:
                    logger.info("analyze {},{}",tt, frame)
                    # stable_evs, np.sum(stable_evs), (stable_th_lev / 100) * window_lev, evs[frame]
                end_ev = frame + 1
                ev = True

            last = (last + 1) % window_lev  # Update pointer to last element

        # Check both
        if ev is True and ess is True:
            end_both = frame + 1
    logger.info("terminate {},{},{},{}",tt,end_ev, end_ess, end_both)
    return end_ev, end_ess, end_both, ess_by_frame, precision_by_frame


def generate_annotation(levs_df, traj, stable_th_lev,window_lev,tt, window_ess, var_th_ess):
    annotations = []
    traj_levs = levs_df.loc[
                    (levs_df["traj_id"] == str(traj.traj_id))
                    & (levs_df["subset"] == traj.subset), "levs"]

    if tt == 0 or tt == 20:
        logger.info("anotation {},{}",tt, traj_levs)
    

    levs = traj_levs.values[0]

    end_ev, end_ess, end_both, ess, p = analyse_trajectory(traj.mdtraj, levs,stable_th_lev,window_lev,tt, window_ess, var_th_ess)
    terminations = [end_ev, end_ess, end_both]
    logger.debug(
        "Trajectory {} stabilized in frame {} (min(LEV, ESS, LEV+ESS))".format(
            traj.traj_id, min(terminations)
        )
    )

    # Save for data frame
    annotations.append(
        [
            traj.traj_id,
            traj.conformation,
            traj.stride,
            traj.n_frames,
            traj.length,
            ess,
            p,
            end_ev,
            end_ess,
            end_both,
        ]
    )

    columns = [
        "traj_id",
        "conformation",
        "stride",
        "n_frames",
        "length",
        "ess",
        "precision",
        "termination_lev",
        "termination_ess",
        "termination_lev_ess",
    ]
    annotation_df = pd.DataFrame(annotations, columns=columns)
    # logger.debug(annotation_df)
    return annotation_df


"""
# Using parallelism, load the trajectories into a dataframe. This function should be able to load multiple trajectories from multiple directories
# and return a single dataframe with all the trajectories. The function should be able to load the trajectories in parallel
# using threads. input_dirs is a list of directories where the trajectories are stored. Each directory contains multiple trajectories.
# the function to load a single trajectory receives the directory and the filename of the trajectory to load.
"""


def load_trajectories(input_dirs):
    trajectory_df = pd.DataFrame()

    with multiprocessing.Pool() as pool:
        # Map the loading of each trajectory to the pool of workers
        results = pool.starmap(
            load_trajectory,
            [
                (directory, filename)
                for directory in input_dirs
                for filename in os.listdir(directory)
            ],
        )
        trajectory_df = pd.concat(results, ignore_index=True)
    return trajectory_df


def generate_annotations(levs_df, trajectory_df, stable_th_lev,window_lev, window_ess, var_th_ess):
    annotations_df = pd.DataFrame()

    # serial for testing
    # results = generate_annotation(levs_df, trajectory_df.iloc[0], stable_th_lev,window_lev,0)
    # return results

    with multiprocessing.Pool() as pool:
        # Map the loading of each trajectory to the pool of workers
        results = pool.starmap(
            generate_annotation,
            [(levs_df, trajectory_df.iloc[tt], stable_th_lev,window_lev,tt, window_ess, var_th_ess) for tt in range(len(trajectory_df))],
        )
        annotations_df = pd.concat(results, ignore_index=True)
    return annotations_df


## -------------------- Main -------------------- ##
if __name__ == "__main__":
    arguments = sys.argv[1:]
    print(sys.argv[0],arguments)
    termination_criterion = arguments[0]
    DEBUG = arguments[1] == "True"

    if not DEBUG:
        output_folder = arguments[7]
        logger.info(f"Saving output to {output_folder}")

    # # Set up logger
    logger.remove(0)
    # logger.add(sys.stderr, level="INFO")
    logger.add(f"{output_folder}/logs.log", level="DEBUG")
    # logger.disable("__main__")

    ## -------------------- Definitions -------------------- ##
    stride = int(arguments[2])  # Stride for loading trajectories
    sampling_rate = float(arguments[3])  # In ps
    to_ns = sampling_rate * stride  # Conversion to ns
    # LEV
    range_min = float(arguments[4])  # LEV range
    range_max = float(arguments[5])  # LEV range
    range_tolerance = float(arguments[6])  # LEV tolerance to remain in range (in %)
    range_th = (range_tolerance / 100) * (range_max - range_min)  # LEV range threshold

    
    # # Load trajectories
    input_dirs = arguments[8].split(',')
    print(input_dirs)
    trajectory_df = load_trajectories(input_dirs)

    # # calculate total of frames in all trajectories
    total_frames = trajectory_df["n_frames"].sum()
    logger.info("Total frames: {}", total_frames)
    logger.debug(trajectory_df)
    
    # # Load eigenvalues
    # # Read the file and create a dataframe
    levs_df = pd.read_pickle(arguments[9])
    

    # # # Generate annotations
    # # Define analysis parameters
    window_levs = [int(x) for x in arguments[10].split(',')]  # Sliding window size for LEV
    stable_th_levs = [int(x) for x in arguments[11].split(',')]  # Minimum number of frames considered in window to decide termination for LEV
    var_th_esss = [int(x) for x in arguments[12].split(',')]  # Minimum average differential variability threshold for ESS (in %)
    window_esss = [int(x) for x in arguments[13].split(',')]  # Sliding window size for ESS
    print(window_levs, stable_th_levs, var_th_esss, window_esss)
    for window_lev in window_levs:
        for stable_th_lev in stable_th_levs:
            for window_ess in window_esss:
                for var_th_ess in var_th_esss:
                    annotations_df = generate_annotations(levs_df, trajectory_df, stable_th_lev, window_lev, window_ess, var_th_ess)
                    sum_lev_sum = annotations_df["termination_lev"].sum()
                    sum_ess_sum = annotations_df["termination_ess"].sum()
                    logger.debug("Sum of LEV terminations: {}", sum_lev_sum)
                    logger.debug("Sum of ESS terminations: {}", sum_ess_sum)
                    # print percentage of frames that terminated
                    percentage_lev = round((sum_lev_sum / total_frames) * 100, 2)
                    percentage_ess = round((sum_ess_sum / total_frames) * 100, 2)
                    if termination_criterion == "termination_lev":
                        file_name = f"{output_folder}/annotations_{termination_criterion}_{percentage_lev}_perc_{stable_th_lev}_window_lev_{window_lev}_window_ess_{window_ess}"
                    else:
                        file_name = f"{output_folder}/annotations_{termination_criterion}_{percentage_ess}_perc_{var_th_ess}_window_lev_{window_lev}_window_ess_{window_ess}"
                    logger.debug("Percentage of frames that terminated by LEV: {}", percentage_lev)
                    logger.debug("Percentage of frames that terminated by ESS: {}", percentage_ess)
                    annotations_df.to_pickle(
                        file_name + ".pkl", protocol=4
                    )  # Protocol 4 necessary for compatibility with Python 3.6 used in MSMBuilder
                    logger.debug("Annotations saved to {}".format(file_name))
        logger.debug(annotations_df)
        logger.info("annotations saved as: {}", file_name)
    
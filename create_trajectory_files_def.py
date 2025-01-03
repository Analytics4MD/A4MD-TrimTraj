import multiprocessing
import sys
import matplotlib.pyplot as plt
import numpy as np
import mdtraj as md
import pyemma.coordinates as coor
from pyemma.util.contexts import settings
from loguru import logger
import pandas as pd

from utils import load_annotation_files, load_trajectories, trim_trajectories
import datetime
import os
import glob

## TODO: move to utils.py file
def save_trajectory(traj, filename):
    traj.save(filename)
    logger.info("Trajectory saved: {}", filename)
    # logger.debug("frames for: {},{}", traj.n_frames, filename)

## TODO: move to utils.py file
def trimming_trajectories(
    trajectories_df,
    annotations_list_df,
    input_dirs,
    stride,
    first_batch,
    trim_first_rand,
    run_last_full,
    termination_criterion,
    trial,
):
    logger.info("Trimming trajectories in parallel")
    with multiprocessing.Pool() as pool:
        # Map the loading of each trajectory to the pool of workers
        results = pool.starmap(
            trim_trajectories,
            [
                (
                    trajectories_df,
                    annotation_df,
                    input_dirs,
                    stride,
                    first_batch,
                    trim_first_rand,
                    run_last_full,
                    termination_criterion,
                    trial,
                )
                for annotation_df in annotations_list_df
            ],
        )
    new_trajectory_file_names = [key for result in results for key in result.keys()]
    results = [item for result in results for item in result.values()]
    # results = [list(result.values()) for result in results]
    logger.debug(
        "Results of trimming names: {},{}",
        new_trajectory_file_names,
        len(new_trajectory_file_names),
    )
    logger.debug("Results of trimming: {} -- {}", (results), len(results))

    # Join the trajectories in the results list and save them to disk in parallel
    with multiprocessing.Pool() as pool:
        trimmed_mdtrajs = pool.map(md.join, results)
    logger.debug(
        "len trimmed: {}, {},{}",
        len(trimmed_mdtrajs),
        type(trimmed_mdtrajs[0]),
        trimmed_mdtrajs,
    )
    logger.debug(
        "frames: {},{}",
        trimmed_mdtrajs[0].n_frames,
        (trimmed_mdtrajs[0].n_frames / 800040),
    )

    if not DEBUG:
        for traj in trimmed_mdtrajs:
            file_name = f"{output_folder}/{termination_criterion}_{round((traj.n_frames*100/5500040),2)}.xtc"
            save_trajectory(traj, file_name)

    return []


def loading_annotation_files(annotation_paths, input_dirs):
    with multiprocessing.Pool() as pool:
        # Map the loading of each trajectory to the pool of workers
        results = pool.starmap(
            load_annotation_files,
            [(annotation_path, input_dirs) for annotation_path in annotation_paths],
        )

    return results


if __name__ == "__main__":
    arguments = sys.argv[1:]
    print(sys.argv[0],arguments)
    termination_criterion = arguments[0]
    DEBUG = arguments[1] == "True"

    if not DEBUG:
        output_folder = arguments[2]
        logger.info(f"Saving output to {output_folder}")

    # # Set up logger
    logger.remove(0)
    # logger.add(sys.stderr, level="INFO")
    logger.add(f"{output_folder}/logs.log", level="DEBUG")
    # logger.disable("__main__")

    ## -*-*-*-*-*-*-*-*-*-*- Loading trajectories -*-*-*-*-*-*-*-*-*-*-
    # # Definitions
    stride = int(arguments[3])  # Stride for loading trajectories
    top_file = arguments[4]  # Topology file for the trajectories
    end = None  # None means load all frames
    # Load full trajectories files
    input_dirs = arguments[5].split(',')
    
    trajectories_df = load_trajectories(input_dirs,stride=stride, end=end)

    # Load full trajectories
    trajs = trajectories_df["mdtraj"].tolist()
    logger.debug("{},{}", trajs, len(trajs))
    logger.info("Trajectories: {}", len(trajs))
    md_trajectories = md.join(trajs)
    logger.info("Join done")

    # # -*-*-*-*-*-*-*-*-*-*- Loading features -*-*-*-*-*-*-*-*-*-*-
    feat = coor.featurizer(top_file)
    feat.add_backbone_torsions()
    # features_ref = coor.load(trajs, feat)
    features_ref = feat.transform(md_trajectories)


    # #  -*-*-*-*-*-*-*-*-*-*- Loading annotations file -*-*-*-*-*-*-*-*-*-*-
    logger.info("Loading annotations")
    annotation_folder = arguments[7]
    annotation_files = os.listdir(annotation_folder)
    annotation_files = glob.glob(os.path.join(annotation_folder, "*.pkl"))
    annotation_files_lev = []
    annotation_files_ess = []
    for file in annotation_files:
        print(file)
        if file.startswith("annotations_termination_lev"):
            annotation_files_lev.append(file)
        else:
            annotation_files_ess.append(file)

    if termination_criterion == "termination_lev":
        annotation_paths = annotation_files_lev
        logger.info("Using lev annotations: {}", annotation_paths)
    else:
        annotation_paths = annotation_files_ess
        logger.info("Using ess annotations: {}", annotation_paths)
    
    annotations_list_df = loading_annotation_files(annotation_paths, input_dirs)
    logger.info("{},{}", annotations_list_df, len(annotations_list_df))
    logger.info("Loaded Annotations: {}", len(annotations_list_df))

    # # # # # -*-*-*-*-*-*-*-*-*-*- Trimming trajectories to termination -*-*-*-*-*-*-*-*-*-*-
    # # # Definitions
    first_batch = int(arguments[8])  # FS1 or FS2 runs first (ie 1 or 2 are trimmed), 0 means trim everything, -1 means dont trim
    trim_first_rand = arguments[9] == "True"  # Terminate at a random point
    trial = arguments[10]  # Identifier of random trial [0, ..., 4]
    run_last_full = arguments[11] == "True"  # Let restarted trajectory run to completion
    logger.info("Trimming trajectories")

    trajectories_filenames = trimming_trajectories(
        trajectories_df,
        annotations_list_df,
        input_dirs,
        stride,
        first_batch,
        trim_first_rand,
        run_last_full,
        termination_criterion,
        trial,
    )

    logger.info("Trajectories trimmed: {}", trajectories_filenames)

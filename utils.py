import os
import mdtraj as md
import pandas as pd
import multiprocessing
import random

from loguru import logger

def load_trajectory(input_dirs, filename,stride=1,end=None):
    sampling_rate = 0.002  # In ps
    to_ns = sampling_rate * stride  # Conversion to ns
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
    if end is not None:
        end = int(end * t.n_frames)
        # end = int((1-end)*t.n_frames)
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
            t if end is None else t[:end],
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

    # print("{} trajectories".format(len(trajectories)))
    return trajectory_df

def load_trajectories(input_dirs,stride=1, end=None):
    trajectory_df = pd.DataFrame()

    with multiprocessing.Pool() as pool:
        # Map the loading of each trajectory to the pool of workers
        results = pool.starmap(
            load_trajectory,
            [
                (directory, filename,stride, end)
                for directory in input_dirs
                for filename in os.listdir(directory)
            ],
        )
        trajectory_df = pd.concat(results, ignore_index=True)
    return trajectory_df

def load_annotation_files(annotation_file, input_dirs):
    dfs = []

    dfs.append(pd.read_pickle(annotation_file))
    subset_1 = [input_dirs[0]] * (len(dfs[0]) // 2)
    subset_2 = [input_dirs[1]] * (len(dfs[0]) // 2)
    subset = subset_1 + subset_2
    dfs[0]["subset"] = subset
    annotation_df = pd.concat(dfs, ignore_index=True)

    # logger.info(annotation_df)

    return annotation_df

def trim_trajectories(trajectories_df, annotation_df, input_dirs,stride, first_batch=0, trim_first_rand=False, run_last_full=False, termination_criterion="termination_lev", trial=0):

    executed_steps = 0
    total_steps = 0
    xyz = []
    random.seed(trial)

    trajectories_df["traj_id"] = trajectories_df["traj_id"].astype(int)

    for tt in range(len(trajectories_df)):
        traj = trajectories_df.iloc[tt]  # sorted_df.iloc[tt]
        total_steps += stride * (traj.n_frames - 1) * 1000

        # if not (traj.subset==input_dirs[0]):
        #    continue

        if (first_batch == 0 or first_batch == 1 and int(traj.traj_id) <= 10) or (
            first_batch == 2 and int(traj.traj_id) > 10
        ):  # Cut to annotation
            if trim_first_rand:
                traj_end = random.randint(0, traj.n_frames - 1)
            else:
                traj_end = annotation_df.loc[
                    (annotation_df["traj_id"] == str(traj.traj_id))
                    & (annotation_df["subset"] == traj.subset),
                    termination_criterion,
                ].values[0]

            logger.debug(
                "Trajectory {} in subset {} cut to {} annotation: now ending {}".format(
                    traj.traj_id, traj.subset, termination_criterion, traj_end
                )
            )
        elif (first_batch == 1 and int(traj.traj_id) > 10) or (
            first_batch == 2 and int(traj.traj_id) <= 10
        ):  # Cut to max steps
            traj_length = traj.n_frames

            if first_batch == 1:
                prev_traj_id = str(int(traj.traj_id) - 10)
            else:
                prev_traj_id = str(int(traj.traj_id) + 10)

            if run_last_full:
                traj_end = traj.n_frames  # To let second trajectory run to completion
            else:
                if trim_first_rand:
                    index_prev = 20 - int(prev_traj_id)
                    print(index_prev)
                    prev_steps = len(xyz[index_prev])
                else:
                    prev_steps = annotation_df.loc[
                        (annotation_df["traj_id"] == prev_traj_id)
                        & (annotation_df["subset"] == traj.subset),
                        termination_criterion,
                    ].values[0]

                traj_end = traj.n_frames - prev_steps

            logger.debug(
                "Trajectory {} in subset {} cut to fit {} steps. {} steps removed from trajectory {}: now ending {}".format(
                    traj.traj_id,
                    traj.subset,
                    traj.n_frames,
                    prev_steps,
                    prev_traj_id,
                    traj_end,
                )
            )

            if traj_end == 0:
                logger.debug("Skipping trajectory")
                continue
        elif first_batch == -1:
            logger.debug(
                "Trajectory {} in subset {} ending at last frame {}".format(
                    traj.traj_id, traj.subset, traj.n_frames
                )
            )
            xyz.append(traj.mdtraj)
            continue
        else:
            logger.debug("Weird case :/")
            continue

        executed_steps += stride * (traj_end - 1) * 1000
        xyz.append(
            traj.mdtraj[:traj_end]
        )  # Copying to xyz for compatibility with example code
    perc = round((executed_steps * 100 / total_steps), 2)
    k = str(perc)+ "_" + termination_criterion
    logger.debug(
        "MD steps executed: {} out of {} ({}%) using {}".format(
            executed_steps, total_steps, perc, termination_criterion
        )
    )
    return {k:xyz}

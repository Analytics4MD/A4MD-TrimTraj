import pyemma
import matplotlib.pyplot as plt
import pandas as pd
import mdtraj as md
import numpy as np
import pyemma
import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import rms



def plot_free_energy(
    cc_x,
    cc_y,
    free_energy_per_cluster,
    free_energy_per_clusterT,
    free_energy_per_clusterT_ess,
    dist_cmap,
    title,
    nstates,
    pcca_sets,
    cols,
    size,
    frames_closest_to_minimum_energy_coor,
    frames_closest_to_minimum_energy_coorT,
    frames_closest_to_minimum_energy_coorT_ess,
):
    """
    Plots the free energy landscape and centroids for different states.

    Parameters:
    cc_x (array-like): X coordinates of the conformational space.
    cc_y (array-like): Y coordinates of the conformational space.
    free_energy_per_cluster (array-like): Free energy values per cluster for the full conformational space.
    free_energy_per_clusterT (array-like): Free energy values per cluster for the trimmed conformational space (LEV).
    free_energy_per_clusterT_ess (array-like): Free energy values per cluster for the trimmed conformational space (ESS).
    dist_cmap (str or colormap): Color map for the contour plot.
    title (str): Title of the plot.
    nstates (int): Number of states.
    pcca_sets (list of arrays): List of arrays containing the indices of the clusters for each state.
    cols (list of str): List of colors for each state.
    size (float): Size of the scatter points.
    frames_closest_to_minimum_energy_coor (array-like): Coordinates of the frames closest to the minimum energy for each state in the full conformational space.
    frames_closest_to_minimum_energy_coorT (array-like): Coordinates of the frames closest to the minimum energy for each state in the trimmed conformational space (LEV).
    frames_closest_to_minimum_energy_coorT_ess (array-like): Coordinates of the frames closest to the minimum energy for each state in the trimmed conformational space (ESS).
    """
    
    fig1, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

    tcf1 = ax1.tricontourf(cc_x, cc_y, free_energy_per_cluster, 200, cmap=dist_cmap)
    fig1.colorbar(tcf1, ax=[ax1, ax2, ax3])
    ax1.set_title("Full")
    for i in range(nstates):
        ax1.scatter(cc_x[pcca_sets[i]], cc_y[pcca_sets[i]], color=cols[i], s=size)
        centroid_label = f"centroid {cols[i]} {i} ({frames_closest_to_minimum_energy_coor[i][0]:.4f}, {frames_closest_to_minimum_energy_coor[i][1]:.2f})"
        ax1.scatter(
            frames_closest_to_minimum_energy_coor[i][0],
            frames_closest_to_minimum_energy_coor[i][1],
            color="white",
            marker="*",
            label=centroid_label,
        )

    tcf2 = ax2.tricontourf(cc_x, cc_y, free_energy_per_clusterT, 200, cmap=dist_cmap)
    ax2.set_title("Trim - LEV")
    for i in range(nstates):
        ax2.scatter(cc_x[pcca_sets[i]], cc_y[pcca_sets[i]], color=cols[i], s=size)
        centroid_label = f"centroid {cols[i]} {i} ({frames_closest_to_minimum_energy_coorT[i][0]:.4f}, {frames_closest_to_minimum_energy_coorT[i][1]:.2f})"
        ax2.scatter(
            frames_closest_to_minimum_energy_coorT[i][0],
            frames_closest_to_minimum_energy_coorT[i][1],
            color="white",
            marker="*",
            label=centroid_label,
        )

    tcf3 = ax3.tricontourf(
        cc_x, cc_y, free_energy_per_clusterT_ess, 200, cmap=dist_cmap
    )
    ax3.set_title("Trim - ESS")
    for i in range(nstates):
        ax3.scatter(cc_x[pcca_sets[i]], cc_y[pcca_sets[i]], color=cols[i], s=size)
        centroid_label = f"centroid {cols[i]} {i} ({frames_closest_to_minimum_energy_coorT_ess[i][0]:.4f}, {frames_closest_to_minimum_energy_coorT_ess[i][1]:.2f})"
        ax3.scatter(
            frames_closest_to_minimum_energy_coorT_ess[i][0],
            frames_closest_to_minimum_energy_coorT_ess[i][1],
            color="white",
            marker="*",
            label=centroid_label,
        )
    fig1.suptitle(title)




def find_frames_closest_to_minimum_energy(
    FES,
    FES_T,
    FES_T_ess,
    M,
    MT,
    MT_ess,
    clust,
    tica_concatenated_full,
    tica_concatenated_trim,
    tica_concatenated_trim_ess,
    nstates,
):
    """
    Finds the frames closest to the minimum energy coordinates for each state.

    Args:
        FES (numpy.ndarray): Free energy surface.
        FES_T (numpy.ndarray): Trimmed free energy surface.
        FES_T_ess (numpy.ndarray): Trimmed essential free energy surface.
        M (object): Metastable object.
        MT (object): Trimmed metastable object.
        MT_ess (object): Trimmed essential metastable object.
        clust (object): Clustering object.
        tica_concatenated_full (numpy.ndarray): Concatenated tICA coordinates.
        tica_concatenated_trim (numpy.ndarray): Trimmed concatenated tICA coordinates.
        tica_concatenated_trim_ess (numpy.ndarray): Trimmed essential concatenated tICA coordinates.
        nstates (int): Number of states.

    Returns:
        tuple: A tuple containing the following:
            - frames_closest_to_minimum_energy_coor (list): List of frames coordinates closest to the minimum energy.
            - frames_closest_to_minimum_energy_coor_T (list): List of frames coordinates closest to the minimum energy (LEV-trimmed).
            - frames_closest_to_minimum_energy_coor_T_ess (list): List of frames coordinates closest to the minimum energy (ESS-trimmed).
            - frames_closest_to_minimum_energy (list): List of frames closest to the minimum energy.
            - frames_closest_to_minimum_energyT (list): List of frames closest to the minimum energy (LEV-trimmed).
            - frames_closest_to_minimum_energyT_ess (list): List of frames closest to the minimum energy (ESS-trimmed).
            - closest_10_frame_index_total (list): List of indexes of the 10 closest frames to the minimum energy coordinates.
            - closest_10_frame_indexT_total (list): List of indexes of the 10 closest frames to the minimum energy coordinates (LEV-trimmed).
            - closest_10_frame_indexT_ess_total (list): List of indexes of the 10 closest frames to the minimum energy coordinates (ESS-trimmed).
            - minimum_energy (list): List of minimum energy values.
            - minimum_energyT (list): List of minimum energy values (LEV-trimmed).
            - minimum_energyT_ess (list): List of minimum energy values (ESS-trimmed).
    """
    pcca_sets = M.metastable_sets
    pcca_setsT = MT.metastable_sets
    pcca_setsT_ess = MT_ess.metastable_sets
    frames_closest_to_minimum_energy_coor = []
    frames_closest_to_minimum_energy_coor_T = []
    frames_closest_to_minimum_energy_coor_T_ess = []
    tICA_first_two = tica_concatenated_full[:, :2]
    frames_closest_to_minimum_energy = []
    tICA_first_twoT = tica_concatenated_trim[:, :2]
    frames_closest_to_minimum_energyT = []
    tICA_first_twoT_ess = tica_concatenated_trim_ess[:, :2]
    frames_closest_to_minimum_energyT_ess = []
    (
        closest_10_frame_index_total,
        closest_10_frame_indexT_total,
        closest_10_frame_indexT_ess_total,
    ) = [], [], []
    minimum_energy, minimum_energyT, minimum_energyT_ess = [], [], []
    comparison_lev=[0,1,4,2,3,5] ## TODO: set as argument from users
    comparison_ess=[2,3,4,1,5,0] ## TODO: set as argument from users

    for state in range(nstates):
        # Gets the indexes of the cluster centers in the state
        indexes_cluster_centers_per_state = pcca_sets[state]
        indexes_cluster_centers_per_stateT = pcca_setsT[state]
        indexes_cluster_centers_per_stateT_ess = pcca_setsT_ess[state]
        # Gets the coordinates of the cluster centers in the state
        coord_cluster_per_state = clust.clustercenters[
            indexes_cluster_centers_per_state
        ]
        coord_cluster_per_stateT = clust.clustercenters[
            indexes_cluster_centers_per_stateT
        ]
        coord_cluster_per_stateT_ess = clust.clustercenters[
            indexes_cluster_centers_per_stateT_ess
        ]

        # Gets the free energy of the cluster centers in the state
        energy_of_cluster_center_per_state = FES[indexes_cluster_centers_per_state]
        energy_of_cluster_center_per_stateT = FES_T[indexes_cluster_centers_per_stateT]
        energy_of_cluster_center_per_stateT_ess = FES_T_ess[
            indexes_cluster_centers_per_stateT_ess
        ]

        # Gets the index of the cluster center with minimum energy in the state
        index_minimum_energy = np.argmin(energy_of_cluster_center_per_state, axis=0)
        minimum_energy.append(np.min(energy_of_cluster_center_per_state))
        index_minimum_energyT = np.argmin(energy_of_cluster_center_per_stateT, axis=0)
        minimum_energyT.append(np.min(energy_of_cluster_center_per_stateT))
        index_minimum_energyT_ess = np.argmin(
            energy_of_cluster_center_per_stateT_ess, axis=0
        )
        minimum_energyT_ess.append(np.min(energy_of_cluster_center_per_stateT_ess))
        ## FINDING THE CLOSEST FRAME TO THE FREE MINIMUM ENERGY COORDINATE
        dist = tICA_first_two - coord_cluster_per_state[index_minimum_energy]
        distances = np.linalg.norm(dist, axis=1)
        closest_frame_index = np.argmin(distances)
        closest_10_frame_index = np.argsort(distances)[:10]
        coor_minimum_energy = tICA_first_two[closest_frame_index]
        frames_closest_to_minimum_energy.append(closest_frame_index)

        dist = tICA_first_twoT - coord_cluster_per_stateT[index_minimum_energyT]
        distances = np.linalg.norm(dist, axis=1)
        closest_frame_indexT = np.argmin(distances)
        closest_10_frame_indexT = np.argsort(distances)[:10]
        coor_minimum_energyT = tICA_first_twoT[closest_frame_indexT]
        frames_closest_to_minimum_energyT.append(closest_frame_indexT)

        dist = (
            tICA_first_twoT_ess
            - coord_cluster_per_stateT_ess[index_minimum_energyT_ess]
        )
        distances = np.linalg.norm(dist, axis=1)
        closest_frame_indexT_ess = np.argmin(distances)
        closest_10_frame_indexT_ess = np.argsort(distances)[:10]
        coor_minimum_energyT_ess = tICA_first_twoT_ess[closest_frame_indexT_ess]
        frames_closest_to_minimum_energyT_ess.append(closest_frame_indexT_ess)

        frames_closest_to_minimum_energy_coor.append(
            # coord_cluster_per_state[index_minimum_energy]
            coor_minimum_energy
        )  # to plot the centers with minimum energy
        frames_closest_to_minimum_energy_coor_T.append(
            # coord_cluster_per_stateT[index_minimum_energyT]
            coor_minimum_energyT
        )  # to plot the centers with minimum energy
        frames_closest_to_minimum_energy_coor_T_ess.append(
            # coord_cluster_per_stateT_ess[index_minimum_energyT_ess]
            coor_minimum_energyT_ess
        )  # to plot the centers with minimum energy
        closest_10_frame_index_total.append(closest_10_frame_index)
        closest_10_frame_indexT_total.append(closest_10_frame_indexT)
        closest_10_frame_indexT_ess_total.append(closest_10_frame_indexT_ess)
        closest_10_frame_index, closest_10_frame_indexT, closest_10_frame_indexT_ess = (
            [],
            [],
            [],
        )

    frames_closest_to_minimum_energy_coor_T = [frames_closest_to_minimum_energy_coor_T[i] for i in comparison_lev]
    frames_closest_to_minimum_energy_coor_T_ess = [frames_closest_to_minimum_energy_coor_T_ess[i] for i in comparison_ess]
    frames_closest_to_minimum_energyT = [frames_closest_to_minimum_energyT[i] for i in comparison_lev]
    frames_closest_to_minimum_energyT_ess = [frames_closest_to_minimum_energyT_ess[i] for i in comparison_ess]
    closest_10_frame_indexT_total = [closest_10_frame_indexT_total[i] for i in comparison_lev]
    closest_10_frame_indexT_ess_total = [closest_10_frame_indexT_ess_total[i] for i in comparison_ess]
    minimum_energyT = [minimum_energyT[i] for i in comparison_lev]
    minimum_energyT_ess = [minimum_energyT_ess[i] for i in comparison_ess]

    return (
        frames_closest_to_minimum_energy_coor,
        frames_closest_to_minimum_energy_coor_T,
        frames_closest_to_minimum_energy_coor_T_ess,
        frames_closest_to_minimum_energy,
        frames_closest_to_minimum_energyT,
        frames_closest_to_minimum_energyT_ess,
        closest_10_frame_index_total,
        closest_10_frame_indexT_total,
        closest_10_frame_indexT_ess_total,
        minimum_energy,
        minimum_energyT,
        minimum_energyT_ess,
    )


def save_frames(folder_path, energy_method, frames_closest, mdtrajectories, prefix):
    """
    Save the frames closest to each state in a given folder.

    Parameters:
    folder_path (str): The path to the folder where the frames will be saved.
    energy_method (str): The energy method used to calculate the frames.
    frames_closest (list): A list of indices representing the frames closest to each state.
    mdtrajectories (list): A list of MDTrajectory objects representing the trajectories.
    prefix (str): The prefix to be added to the saved file names.

    Returns:
    frames_minimum (list): A list of MDTrajectory objects representing the saved frames.
    frames_minimum_files (list): A list of file paths representing the saved frames.
    """
    frames_minimum = []
    frames_minimum_files = []

    for state, frame in enumerate(frames_closest):
        frame_in_trajectory = mdtrajectories[frame]
        file_name = os.path.join(
            folder_path, f"{prefix}_state_{state}_frame_{frame}_{energy_method}.pdb"
        )
        frames_minimum_files.append(file_name)
        frames_minimum.append(frame_in_trajectory)
        frame_in_trajectory.save(file_name)

    return frames_minimum, frames_minimum_files


def save_10_frames(folder_path, energy_method, frames_closest, mdtrajectories, prefix):
    """
    Save the 10 frames closest to each state in the frames_closest list.

    Parameters:
    - folder_path (str): The path to the folder where the frames will be saved.
    - energy_method (str): The energy method used to calculate the frames.
    - frames_closest (list): A list of lists containing the indices of the 10 closest frames for each state.
    - mdtrajectories (list): A list of mdtraj.Trajectory objects representing the MD trajectories.
    - prefix (str): The prefix to be added to the file names.

    Returns:
    - frames_10_minimum_total (list): A list of lists containing the 10 frames closest to each state.
    - frames_10_minimum_files_total (list): A list of lists containing the file names of the saved frames.
    """
    frames_10_minimum = []
    frames_10_minimum_files = []
    frames_10_minimum_total = []
    frames_10_minimum_files_total = []

    for state, frames_10 in enumerate(frames_closest):
        for i, frame in enumerate(frames_10):
            frame_in_trajectory = mdtrajectories[frame]
            file_name = os.path.join(
                folder_path,
                f"{prefix}_{energy_method}_state_{state}_{i}_frame_{frame}.pdb",
            )
            frames_10_minimum_files.append(file_name)
            frames_10_minimum.append(frame_in_trajectory)
            frame_in_trajectory.save(file_name)
        frames_10_minimum_total.append(frames_10_minimum)
        frames_10_minimum_files_total.append(frames_10_minimum_files)
        frames_10_minimum = []
        frames_10_minimum_files = []

    return frames_10_minimum_total, frames_10_minimum_files_total


def create_pyemma_energy_plots(
    tica_concatenated_full,
    tica_concatenated_trim,
    tica_concatenated_trim_ess,
    M,
    MT,
    MT_ess,
    pcca_sets,
    frames_closest_to_minimum_energy_coor,
    frames_closest_to_minimum_energy_coor_T,
    frames_closest_to_minimum_energy_coor_T_ess,
    nstates,
    cc_x,
    cc_y,
    cols,
    size,
    title,
):
    """
    Create PyEMMA energy plots.

    Parameters:
    - tica_concatenated_full: numpy.ndarray
        TICA concatenated full array.
    - tica_concatenated_trim: numpy.ndarray
        TICA concatenated LEV-trimmed array.
    - tica_concatenated_trim_ess: numpy.ndarray
        TICA concatenated ESS-trimmed array.
    - M: pyemma.msm.markov_model
        MSM object for full data.
    - MT: pyemma.msm.markov_model
        MSM object for LEV-trimmmed data.
    - MT_ess: pyemma.msm.markov_model
        MSM object for ESS-trimmed data.
    - pcca_sets: list of numpy.ndarray
        List of PCCA sets.
    - frames_closest_to_minimum_energy_coor: list of numpy.ndarray
        List of frames closest to minimum energy coordinates for each centroid in the full data.
    - frames_closest_to_minimum_energy_coor_T: list of numpy.ndarray
        List of frames closest to minimum energy coordinates for each centroid in the trim data.
    - frames_closest_to_minimum_energy_coor_T_ess: list of numpy.ndarray
        List of frames closest to minimum energy coordinates for each centroid in the trim ESS data.
    - nstates: int
        Number of states.
    - cc_x: numpy.ndarray
        X coordinates of the centroids.
    - cc_y: numpy.ndarray
        Y coordinates of the centroids.
    - cols: list of str
        List of colors for each centroid.
    - size: int
        Size of the scatter points.
    - title: str
        Title of the figure.
    """
    fig, axs = plt.subplots(figsize=(20, 8), sharey=True, ncols=3, sharex=True)
    (ax1, ax2, ax3) = axs

    misc_energy = pyemma.plots.plot_free_energy(
        tica_concatenated_full[:, 0],
        tica_concatenated_full[:, 1],
        # weights=np.concatenate(M.trajectory_weights()),
        cmap="viridis",
        nbins=200,
        cbar=True,
        cbar_label="Full - Free energy/kT",
        zorder=0,
        alpha=1,
        ax=ax1,
        vmin=0,
        vmax=10,
    )
    ax1.grid()
    for i in range(nstates):
        ax1.scatter(cc_x[pcca_sets[i]], cc_y[pcca_sets[i]], color=cols[i], s=size)
        centroid_label = f"centroid {cols[i]} {i} ({frames_closest_to_minimum_energy_coor[i][0]:.4f}, {frames_closest_to_minimum_energy_coor[i][1]:.2f})"
        ax1.scatter(
            frames_closest_to_minimum_energy_coor[i][0],
            frames_closest_to_minimum_energy_coor[i][1],
            color="white",
            marker="*",
            label=centroid_label,
        )
        ax1.legend()

    misc_energyT = pyemma.plots.plot_free_energy(
        tica_concatenated_trim[:, 0],
        tica_concatenated_trim[:, 1],
        # weights=np.concatenate(MT.trajectory_weights()),
        cmap="viridis",
        nbins=200,
        cbar=True,
        cbar_label="Trim - Free energy/kT",
        zorder=0,
        alpha=1,
        ax=ax2,
        vmin=0,
        vmax=10,
    )
    ax2.grid()
    for i in range(nstates):
        ax2.scatter(cc_x[pcca_sets[i]], cc_y[pcca_sets[i]], color=cols[i], s=size)
        centroid_label = f"centroid {i} ({frames_closest_to_minimum_energy_coor_T[i][0]:.4f}, {frames_closest_to_minimum_energy_coor_T[i][1]:.2f})"
        ax2.scatter(
            frames_closest_to_minimum_energy_coor_T[i][0],
            frames_closest_to_minimum_energy_coor_T[i][1],
            color="white",
            marker="*",
            label=centroid_label,
        )
        ax2.legend()

    misc_energyT_ess = pyemma.plots.plot_free_energy(
        tica_concatenated_trim_ess[:, 0],
        tica_concatenated_trim_ess[:, 1],
        # weights=np.concatenate(MT_ess.trajectory_weights()),
        cmap="viridis",
        nbins=200,
        cbar=True,
        cbar_label="Trim ESS - Free energy/kT",
        zorder=0,
        alpha=1,
        ax=ax3,
        vmin=0,
        vmax=10,
    )
    ax3.grid()
    for i in range(nstates):
        ax3.scatter(cc_x[pcca_sets[i]], cc_y[pcca_sets[i]], color=cols[i], s=size)
        centroid_label = f"centroid {i} ({frames_closest_to_minimum_energy_coor_T_ess[i][0]:.4f}, {frames_closest_to_minimum_energy_coor_T_ess[i][1]:.2f})"
        ax3.scatter(
            frames_closest_to_minimum_energy_coor_T_ess[i][0],
            frames_closest_to_minimum_energy_coor_T_ess[i][1],
            color="white",
            marker="*",
            label=centroid_label,
        )
        ax3.legend()
    fig.suptitle(title)


def computeRMSD(u1, u2, nstates, selection):
    """
    Compute the root mean square deviation (RMSD) between two sets of coordinates.

    Parameters:
    u1 (Universe): The first set of coordinates.
    u2 (Universe): The second set of coordinates.
    nstates (int): The number of states.
    selection (str): The selection string for atoms to include in the calculation.

    Returns:
    float: The RMSD value.
    """
    return np.sqrt(
        (
            (
                (
                    u1.select_atoms(selection).positions
                    - u2.select_atoms(selection).positions
                )
                ** 2
            )
            * 3
        ).mean()
    )


def calculate_rmsd_values_s(
    full_u,
    full_p,
    full_s,
    trim_u,
    trim_p,
    trim_s,
    full_u_f,
    full_p_f,
    full_s_f,
    trim_u_f,
    trim_p_f,
    trim_s_f,
    labels,
    method,
    nstates,
    selection,
    title,
    comparison_list,
):
    """
    Calculate and plot RMSD values for different states.

    Args:
        full_u (list): List of umbrella sampling values in the full trajectories for each state.
        full_p (list): List of Pyemma values in the full trajectories values for each state.
        full_s (list): List of Boltzmann equation values in the full trajectories for each state.
        trim_u (list): List of umbrella sampling values in the trimmed trajectories for each state.
        trim_p (list): List of Pyemma values in the trimmed trajectories for each state.
        trim_s (list): List of Boltzmann equation values in the trimmed trajectories for each state.
        full_u_f (list): List of umbrella sampling files in the full trajectories
        full_p_f (list): List of Pyemma files in the full trajectories files for each state.
        full_s_f (list): List of Boltzmann equation files in the full trajectories for each state.
        trim_u_f (list): List of umbrella sampling files in the trimmed trajectories for each state.
        trim_p_f (list): List of Pyemma files in the trimmed trajectories for each state.
        trim_s_f (list): List of Boltzmann equation files in the trimmed trajectories for each state.
        labels (list): List of labels for the heatmap.
        method (str): Method for calculating RMSD ("mdtraj", "compute", or "mdanalysis").
        nstates (int): Number of states.
        selection (str): Atom selection for RMSD calculation.
        title (str): Title for the plot.
        comparison_list (list): List of indices for comparison.

    Returns:
        None
    """
    fig, axs = plt.subplots(
        figsize=(16, 10), sharey=True, ncols=3, nrows=2, sharex=True
    )
    total_values = []

    for state, state_trim in zip(range(nstates), comparison_list):
        frames_in_state = [
            full_u[state],
            full_p[state],
            full_s[state],
            trim_u[state_trim],
            trim_p[state_trim],
            trim_s[state_trim],
        ]
        frames_in_state_f = [
            full_u_f[state],
            full_p_f[state],
            full_s_f[state],
            trim_u_f[state_trim],
            trim_p_f[state_trim],
            trim_s_f[state_trim],
        ]

        rmsd_values = np.zeros((len(frames_in_state), len(frames_in_state)))

        for i in range(len(frames_in_state)):
            for j in range(len(frames_in_state)):
                if method == "mdtraj":
                    rmsd_values[i, j] = md.rmsd(frames_in_state[i], frames_in_state[j])
                elif method == "compute":
                    u1 = mda.Universe(frames_in_state_f[i])
                    u2 = mda.Universe(frames_in_state_f[j])
                    rmsd_values[i, j] = computeRMSD(u1, u2, nstates, selection)
                elif method == "mdanalysis":
                    u1 = mda.Universe(frames_in_state_f[i])
                    u2 = mda.Universe(frames_in_state_f[j])
                    v = rms.rmsd(
                        u1.select_atoms(selection).positions,
                        u2.select_atoms(selection).positions,
                        superposition=True,
                    )
                    rmsd_values[i, j] = v
        total_values.append(rmsd_values)

        vmin = np.min(total_values)
        vmax = np.max(total_values)
    for state in range(nstates):
        df = pd.DataFrame(total_values[state])
        ax = axs[state // 3, state % 3]

        # Create the heatmap using seaborn
        im = sns.heatmap(
            df,
            cmap="viridis",
            ax=ax,
            cbar=True,
            annot=True,
            fmt=".2f",
            cbar_kws={"label": "RMSD"},
            vmin=vmin,
            vmax=vmax,
        )

        # Add title to the heatmap
        ax.set_title(f"State {state}")

        # Set y and x labels
        ax.set_yticklabels(labels)
        ax.set_xticklabels(labels)
    fig.suptitle(title)
    plt.show()


def calculate_average_rmsd(frames_list, frames_files, method, nstates, selection):
    """
    Calculate the average RMSD (Root Mean Square Deviation) for each state in a list of frames.

    Parameters:
    - frames_list (list): A list of frames for each state.
    - frames_files (list): A list of file paths for each frame.
    - method (str): The method to use for calculating RMSD ("mdtraj", "compute", or "mdanalysis").
    - nstates (int): The number of states.
    - selection (str): The selection of atoms to consider for RMSD calculation.

    Returns:
    - rmsd_list (list): A list of average RMSD values for each state.
    """
    rmsd_list = []
    for state, frames in enumerate(frames_list):
        rmsd_sum = 0
        for frame in range(1, len(frames)):
            reference_frame = frames[0]
            reference_frame_file = frames_files[state][0]
            if method == "mdtraj":
                frames[frame].center_coordinates()
                reference_frame.center_coordinates()
                rmsd = md.rmsd(frames[frame], reference_frame, precentered=True).item()
            elif method == "compute":
                u1 = mda.Universe(frames_files[state][frame])
                u2 = mda.Universe(reference_frame_file)
                rmsd = computeRMSD(u1, u2, nstates, selection)
            elif method == "mdanalysis":
                u1 = mda.Universe(frames_files[state][frame])
                u2 = mda.Universe(reference_frame_file)
                rmsd = rms.rmsd(
                    u1.select_atoms(selection).positions,
                    u2.select_atoms(selection).positions,
                    superposition=True,
                    center=True,
                )
            rmsd_sum += rmsd
        average_rmsd = rmsd_sum / (len(frames) - 1)
        rmsd_list.append(average_rmsd)
    return rmsd_list

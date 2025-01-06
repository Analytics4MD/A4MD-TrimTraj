import pyemma
import matplotlib.pyplot as plt
import pandas as pd
import pyemma.coordinates as coor
import mdtraj as md
import numpy as np
import matplotlib as mpl
import pyemma
import matplotlib.pyplot as plt
import pandas as pd
import os
import seaborn as sns
import MDAnalysis as mda
from MDAnalysis.analysis import rms


def plot_free_energy(
    cc_x=None,
    cc_y=None,
    cc_x_lev=None,
    cc_y_lev=None,
    cc_x_ess=None,
    cc_y_ess=None,
    free_energy_per_cluster=None,
    free_energy_per_clusterT=None,
    free_energy_per_clusterT_ess=None,
    dist_cmap=None,
    title=None,
    nstates=None,
    pcca_sets=None,
    pcca_setsT=None,
    pcca_setsT_ess=None,
    cols=None,
    size=None,
    frames_closest_to_minimum_energy_coor=None,
    frames_closest_to_minimum_energy_coorT=None,
    frames_closest_to_minimum_energy_coorT_ess=None,
    no_points=False,
    full_title="Full",
    lev_title="Trim - LEV",
    ess_title="Trim - ESS",
):
    if pcca_setsT is None:
        pcca_setsT = pcca_sets
    if pcca_setsT_ess is None:
        pcca_setsT_ess = pcca_sets
    fig1, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 5))
    vmax = np.max(free_energy_per_cluster)
    tcf1 = ax1.tricontourf(
        cc_x,
        cc_y,
        free_energy_per_cluster,
        200,
        cmap=dist_cmap,
        alpha=0.9,
        linestyles="dotted",
        linewidths=1,
        antialiased=False,  # , vmin=0, vmax=vmax
    )
    cbar_ = fig1.colorbar(tcf1, ax=[ax1, ax2, ax3])
    cbar_.set_label("Free energy / kT", fontsize=14)
    ax1.set_title(full_title, fontsize=20)
    ax1.set_xlabel("TIC 1", fontsize=14)
    ax1.tick_params(
        axis="both", which="major", labelsize=14
    )  # Set ticks fontsize to 14

    ax1.set_ylabel("TIC 2", fontsize=14)
    if not no_points:
        for i in range(nstates):
            ax1.scatter(cc_x[pcca_sets[i]], cc_y[pcca_sets[i]], color=cols[i], s=size)
            centroid_label = f"centroid {cols[i]} {i} ({frames_closest_to_minimum_energy_coor[i][0]:.4f}, {frames_closest_to_minimum_energy_coor[i][1]:.2f})"
            ax1.scatter(
                frames_closest_to_minimum_energy_coor[i][0],
                frames_closest_to_minimum_energy_coor[i][1],
                color="white",
                # color=cols[i],
                marker="*",
                label=centroid_label,
            )
    ax1.grid()
    ax1.set_xlim(-1, 1.8)  # Set x limits

    # ax1.legend()
    if cc_x_lev is None:
        cc_x_lev = cc_x
        cc_y_lev = cc_y
    tcf2 = ax2.tricontourf(
        cc_x_lev,
        cc_y_lev,
        free_energy_per_clusterT,
        200,
        cmap=dist_cmap,
        alpha=0.9,
        linestyles="dotted",
        linewidths=1,
        antialiased=False,  # , vmin=0, vmax=vmax
    )
    ax2.set_title(lev_title, fontsize=20)
    ax2.set_xlabel("TIC 1", fontsize=14)
    ax2.tick_params(
        axis="both", which="major", labelsize=14
    )  # Set ticks fontsize to 14

    # ax2.set_ylabel("TIC 2", fontsize= 14)
    if not no_points:
        for i in range(nstates):
            ax2.scatter(
                cc_x_lev[pcca_setsT[i]], cc_y_lev[pcca_setsT[i]], color=cols[i], s=size
            )
            centroid_label = f"centroid {cols[i]} {i} ({frames_closest_to_minimum_energy_coorT[i][0]:.4f}, {frames_closest_to_minimum_energy_coorT[i][1]:.2f})"
            ax2.scatter(
                frames_closest_to_minimum_energy_coorT[i][0],
                frames_closest_to_minimum_energy_coorT[i][1],
                color="white",
                # color=cols[i],
                marker="*",
                label=centroid_label,
            )
    ax2.grid()
    ax2.set_xlim(-1, 1.8)  # Set x limits

    # ax2.legend()
    if cc_x_ess is None:
        cc_x_ess = cc_x
        cc_y_ess = cc_y
    tcf3 = ax3.tricontourf(
        cc_x_ess,
        cc_y_ess,
        free_energy_per_clusterT_ess,
        200,
        cmap=dist_cmap,
        alpha=0.9,
        linestyles="dotted",
        linewidths=1,
        antialiased=False,  # , vmin=0, vmax=vmax
    )
    ax3.set_title(ess_title, fontsize=20)
    ax3.set_xlabel("TIC 1", fontsize=14)
    ax3.tick_params(
        axis="both", which="major", labelsize=14
    )  # Set ticks fontsize to 14

    ax3.set_ylabel("TIC 2", fontsize=14)
    if not no_points:
        for i in range(nstates):
            ax3.scatter(
                cc_x_ess[pcca_setsT_ess[i]],
                cc_y_ess[pcca_setsT_ess[i]],
                color=cols[i],
                s=size,
            )
            centroid_label = f"centroid {cols[i]} {i} ({frames_closest_to_minimum_energy_coorT_ess[i][0]:.4f}, {frames_closest_to_minimum_energy_coorT_ess[i][1]:.2f})"
            ax3.scatter(
                frames_closest_to_minimum_energy_coorT_ess[i][0],
                frames_closest_to_minimum_energy_coorT_ess[i][1],
                color="white",
                # color=cols[i],
                marker="*",
                label=centroid_label,
            )
    ax3.grid()
    ax3.set_xlim(-1, 1.8)  # Set x limits

    # ax3.legend()
    fig1.suptitle(title)


def find_frames_closest_to_minimum_energy(
    FES=None,
    FES_T=None,
    FES_T_ess=None,
    M=None,
    MT=None,
    MT_ess=None,
    clust=None,
    clustT=None,
    clustT_ess=None,
    tica_concatenated_full=None,
    tica_concatenated_trim_lev=None,
    tica_concatenated_trim_ess=None,
    nstates=None,
):
    pcca_sets = M.metastable_sets
    pcca_setsT = MT.metastable_sets
    pcca_setsT_ess = MT_ess.metastable_sets
    frames_closest_to_minimum_energy_coor = []
    frames_closest_to_minimum_energy_coor_T = []
    frames_closest_to_minimum_energy_coor_T_ess = []
    tICA_first_two = tica_concatenated_full[:, :2]
    frames_closest_to_minimum_energy = []
    tICA_first_twoT = tica_concatenated_trim_lev[:, :2]
    frames_closest_to_minimum_energyT = []
    tICA_first_twoT_ess = tica_concatenated_trim_ess[:, :2]
    frames_closest_to_minimum_energyT_ess = []
    (
        closest_10_frame_index_total,
        closest_10_frame_indexT_total,
        closest_10_frame_indexT_ess_total,
    ) = [], [], []
    minimum_energy, minimum_energyT, minimum_energyT_ess = [], [], []
    if clustT is None:
        clustT = clust
        clustT_ess = clust

    for state in range(nstates):
        # Gets the indexes of the cluster centers in the state
        indexes_cluster_centers_per_state = pcca_sets[state]
        # print("pccasets",pcca_sets[state])
        ## TODO: fix getting the pcca by automation
        indexes_cluster_centers_per_stateT = pcca_sets[state]
        indexes_cluster_centers_per_stateT_ess = pcca_sets[state]
        # Gets the coordinates of the cluster centers in the state
        coord_cluster_per_state = clust.clustercenters[
            indexes_cluster_centers_per_state
        ]
        coord_cluster_per_stateT = clustT.clustercenters[
            indexes_cluster_centers_per_stateT
        ]
        coord_cluster_per_stateT_ess = clustT_ess.clustercenters[
            indexes_cluster_centers_per_stateT_ess
        ]
        energy_of_cluster_center_per_state = FES[indexes_cluster_centers_per_state]
        # print("indexes", indexes_cluster_centers_per_stateT)
        energy_of_cluster_center_per_stateT = FES_T[indexes_cluster_centers_per_stateT]
        energy_of_cluster_center_per_stateT_ess = FES_T_ess[
            indexes_cluster_centers_per_stateT_ess
        ]
        if np.any(energy_of_cluster_center_per_stateT == 0):
            print("lev", state)
        
        # Gets the index of the cluster center with minimum energy in the state
        index_minimum_energy = np.argmin(energy_of_cluster_center_per_state, axis=0)
        minimum_energy.append(np.min(energy_of_cluster_center_per_state))
        # print("state", state, "energy_per_cluster_center", energy_of_cluster_center_per_stateT)
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
    fig, axs = plt.subplots(figsize=(20, 8), sharey=True, ncols=3, sharex=True)
    (ax1, ax2, ax3) = axs

    misc_energy = pyemma.plots.plot_free_energy(
        tica_concatenated_full[:, 0],
        tica_concatenated_full[:, 1],
        cmap="nipy_spectral",
        nbins=200,
        cbar=False,
        cbar_label="Full - Free energy/kT",
        zorder=0,
        alpha=1,
        ax=ax1,
        vmin=0,
        vmax=10,
        legacy=False,
    )
    cbar_ = fig.colorbar(misc_energy[2]["mappable"], ax=[ax1, ax2, ax3])
    cbar_.set_label("Free energy / kT", fontsize=14)
    ax1.grid()

    for i in range(nstates):
        # ax1.scatter(cc_x[pcca_sets[i]], cc_y[pcca_sets[i]], color=cols[i], s=size)
        centroid_label = f"centroid {cols[i]} {i} ({frames_closest_to_minimum_energy_coor[i][0]:.4f}, {frames_closest_to_minimum_energy_coor[i][1]:.2f})"
        ax1.scatter(
            frames_closest_to_minimum_energy_coor[i][0],
            frames_closest_to_minimum_energy_coor[i][1],
            # color="white",
            color=cols[i],
            marker="*",
            label=centroid_label,
        )
        ax1.legend()

    misc_energyT = pyemma.plots.plot_free_energy(
        tica_concatenated_trim[:, 0],
        tica_concatenated_trim[:, 1],
        cmap="nipy_spectral",
        nbins=200,
        cbar=False,
        cbar_label="Trim - Free energy/kT",
        zorder=0,
        alpha=1,
        ax=ax2,
        vmin=0,
        vmax=10,
    )
    ax2.grid()
    for i in range(nstates):
        # ax2.scatter(cc_x[pcca_sets[i]], cc_y[pcca_sets[i]], color=cols[i], s=size)
        centroid_label = f"centroid {i} ({frames_closest_to_minimum_energy_coor_T[i][0]:.4f}, {frames_closest_to_minimum_energy_coor_T[i][1]:.2f})"
        ax2.scatter(
            frames_closest_to_minimum_energy_coor_T[i][0],
            frames_closest_to_minimum_energy_coor_T[i][1],
            # color="white",
            color=cols[i],
            marker="*",
            label=centroid_label,
        )
        ax2.legend()

    misc_energyT_ess = pyemma.plots.plot_free_energy(
        tica_concatenated_trim_ess[:, 0],
        tica_concatenated_trim_ess[:, 1],
        # weights=np.concatenate(MT_ess.trajectory_weights()),
        cmap="nipy_spectral",
        nbins=200,
        cbar=False,
        cbar_label="Trim ESS - Free energy/kT",
        zorder=0,
        alpha=1,
        ax=ax3,
        vmin=0,
        vmax=10,
    )
    ax3.grid()
    for i in range(nstates):
        # ax3.scatter(cc_x[pcca_sets[i]], cc_y[pcca_sets[i]], color=cols[i], s=size)
        centroid_label = f"centroid {i} ({frames_closest_to_minimum_energy_coor_T_ess[i][0]:.4f}, {frames_closest_to_minimum_energy_coor_T_ess[i][1]:.2f})"
        ax3.scatter(
            frames_closest_to_minimum_energy_coor_T_ess[i][0],
            frames_closest_to_minimum_energy_coor_T_ess[i][1],
            # color="white",
            color=cols[i],
            marker="*",
            label=centroid_label,
        )
        ax3.legend()
    fig.suptitle(title)



def computeRMSD(u1, u2, nstates, selection):
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
    full_u=None,
    full_p=None,
    full_s=None,
    trim_u=None,
    trim_p=None,
    trim_s=None,
    full_u_f=None,
    full_p_f=None,
    full_s_f=None,
    trim_u_f=None,
    trim_p_f=None,
    trim_s_f=None,
    labels=None,
    method=None,
    nstates=None,
    selection=None,
    title=None,
    comparison_list=None,
):
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
        ax.set_title(f"State {state}", fontsize=20)

        # Set y and x labels
        ax.set_yticklabels(labels)
        ax.set_xticklabels(labels)
    fig.suptitle(title)
    plt.show()


def calculate_average_rmsd(frames_list, frames_files, method, nstates, selection):
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


def get_histogram(xall, yall, nbins=100, weights=None, avoid_zero_count=False):
    z, xedge, yedge = np.histogram2d(xall, yall, bins=nbins, weights=weights)
    x = 0.5 * (xedge[:-1] + xedge[1:])
    y = 0.5 * (yedge[:-1] + yedge[1:])
    if avoid_zero_count:
        z = np.maximum(z, np.min(z[z.nonzero()]))
    return x, y, z.T  # transpose to match x/y-directions


def _to_density(z):
    """Normalize histogram counts."""
    return z / float(z.sum())


def plot_density(
    xall,
    yall,
    ax=None,
    cmap=None,
    ncontours=100,
    vmin=None,
    vmax=None,
    levels=None,
    cbar=True,
    cax=None,
    cbar_label="Density",
    cbar_orientation="vertical",
    logscale=False,
    nbins=100,
    weights=None,
    avoid_zero_count=False,
    min_v=None,
    max_v=None,
    **kwargs,
):
    """Plot a two-dimensional density map using a histogram of
    scattered data.
    """
    x, y, z = get_histogram(
        xall, yall, nbins=nbins, weights=weights, avoid_zero_count=avoid_zero_count
    )
    pi = _to_density(z)
    # print(pi.min(), pi.max())
    pi = np.ma.masked_where(pi <= 0, pi)

    if logscale:
        from matplotlib.colors import LogNorm

        # print(vmin, vmax)
        norm = LogNorm(vmin=vmin, vmax=vmax)
        if levels is None:
            levels = np.logspace(
                np.floor(np.log10(pi.min())), np.ceil(np.log10(pi.max())), ncontours + 1
            )
            values = np.logspace(
                np.floor(np.log10(pi.min())), np.ceil(np.log10(pi.max()))
            )
            # minmax = (np.floor(np.log10(pi.min())), np.ceil(np.log10(pi.max())))
            min_v = min(values)
            max_v = max(values)
            # print(min_v, max_v)
    else:
        norm = None
    fig, ax, misc = plot_map(
        x,
        y,
        pi,
        ax=ax,
        cmap=cmap,
        ncontours=ncontours,
        vmin=vmin,
        vmax=vmax,
        levels=levels,
        cbar=cbar,
        cax=cax,
        cbar_label=cbar_label,
        cbar_orientation=cbar_orientation,
        norm=norm,
        **kwargs,
    )
    if cbar and logscale:
        from matplotlib.ticker import LogLocator

        misc["cbar"].set_ticks(LogLocator(base=10.0, subs=range(10)))
    # print(np.nonzero(pi)[0])
    return fig, ax, misc, min_v, max_v, levels


def _prune_kwargs(kwargs):
    """Remove non-allowed keys from a kwargs dictionary."""
    allowed_keys = [
        "corner_mask",
        "alpha",
        "locator",
        "extend",
        "xunits",
        "yunits",
        "antialiased",
        "nchunk",
        "hatches",
        "zorder",
    ]
    ignored = [key for key in kwargs.keys() if key not in allowed_keys]
    for key in ignored:
        print(
            "{}={} is not an allowed optional parameter and will" " be ignored".format(
                key, kwargs[key]
            )
        )
        kwargs.pop(key, None)
    return kwargs


def plot_free_energy_pyemma(
    xall,
    yall,
    weights=None,
    ax=None,
    nbins=100,
    ncontours=100,
    offset=-1,
    avoid_zero_count=False,
    minener_zero=True,
    kT=1.0,
    vmin=None,
    vmax=None,
    cmap="nipy_spectral",
    cbar=True,
    cbar_label="free energy / kT",
    cax=None,
    levels=None,
    legacy=True,
    ncountours=None,
    cbar_orientation="vertical",
    **kwargs,
):
    x, y, z = get_histogram(
        xall, yall, nbins=nbins, weights=weights, avoid_zero_count=avoid_zero_count
    )
    f = _to_free_energy(z, minener_zero=minener_zero) * kT
    nonzero = f.nonzero()
    fig, ax, misc = plot_map(
        x,
        y,
        f,
        ax=ax,
        cmap=cmap,
        ncontours=ncontours,
        vmin=vmin,
        vmax=vmax,
        levels=levels,
        cbar=cbar,
        cax=cax,
        cbar_label=cbar_label,
        cbar_orientation=cbar_orientation,
        norm=None,
        **kwargs,
    )
    if legacy:
        return fig, ax
    return fig, ax, misc
    return f


def _to_free_energy(z, minener_zero=False):
    pi = _to_density(z)
    free_energy = np.inf * np.ones(shape=z.shape)
    nonzero = pi.nonzero()
    free_energy[nonzero] = -np.log(pi[nonzero])
    if minener_zero:
        free_energy[nonzero] -= np.min(free_energy[nonzero])
    return free_energy


def plot_map(
    x,
    y,
    z,
    ax=None,
    cmap=None,
    ncontours=100,
    vmin=None,
    vmax=None,
    levels=None,
    cbar=True,
    cax=None,
    cbar_label=None,
    cbar_orientation="vertical",
    norm=None,
    **kwargs,
):
    """Plot a two-dimensional map from data on a grid."""
    import matplotlib.pyplot as _plt

    if ax is None:
        fig, ax = _plt.subplots()
    else:
        fig = ax.get_figure()
    mappable = ax.contourf(
        x,
        y,
        z,
        ncontours,
        norm=norm,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        levels=levels,
        **_prune_kwargs(kwargs),
    )
    misc = dict(mappable=mappable)
    if cbar_orientation not in ("horizontal", "vertical"):
        raise ValueError('cbar_orientation must be "horizontal" or "vertical"')
    if cbar:
        if cax is None:
            cbar_ = fig.colorbar(mappable, ax=ax, orientation=cbar_orientation, pad=0.1)
        else:
            cbar_ = fig.colorbar(
                mappable, cax=cax, orientation=cbar_orientation, pad=0.1
            )
        if cbar_label is not None:
            cbar_.set_label(cbar_label)
        misc.update(cbar=cbar_)
    ax.set_xlabel("TIC 1", fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=14)  # Set ticks fontsize to 14

    ax.set_ylabel("TIC 2", fontsize=14)
    return fig, ax, misc

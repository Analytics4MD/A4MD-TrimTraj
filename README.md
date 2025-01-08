<h1 align="center">  
  A4MD - TrimTraj
  <h4 align="center">

  <a href="https://analytics4md.org/"><img src="https://avatars.githubusercontent.com/u/32650548?s=200&v=4"/></a>

  </h4>
</h1>

<p align="center">
  <a href="#about">About</a> •
  <a href="#prerequisites">Prerequisites</a> •
  <a href="#dependencies">Dependencies</a> •
  <a href="#installation">Installation</a> •
  <a href="#generating-early-terminated-trajectories">Generating early-terminated trajectories</a> •
  <a href="#related-publications">Publications</a> •
  <a href="#copyright-and-license">Copyright and License</a>
</p>

## Authors
Jack Marquez, Silvina Caino-Lores, Michel Cuendet, Trilce Estrada, Ewa Deelman, Harel Weinstein, and Michela Taufer.

## About A4MD

The project's harnessed knowledge of molecular structures' transformations at runtime can be used to steer simulations to more promising areas of the simulation space, identify the data that should be written to congested parallel file systems, and index generated data for retrieval and post-simulation analysis. Supported by this knowledge, molecular dynamics workflows such as replica exchange simulations, Markov state models, and the string method with swarms of trajectories can be executed from the outside (i.e., without reengineering the molecular dynamics code) 

## About A4MD - Validation of conformational space

We validate the A4MD framework capability for early termination and assess whether A4MD early trimmed MD simulations cover the conformational space as effectively as the full simulation.

---
## Prerequisites

In order to use this Jupyter Notebook, your system should have the following installed:
- Anaconda
- python 3.8.18

**NOTE:** It is important to have the simulation data to work with the Jupyter Notebook. The entire dataset with all the trajectories can be found [here](https://doi.org/10.7910/DVN/ML5607)

---
## Dependencies

The framework is also built on top the following third-party libraries: 
- PyEmma
- MDtraj
- MDAnalysis

---
## Installation

Here is the extensive installation instructions. Please make sure the all the prerequisites are satisfied before proceeding the following steps.
Make sure you are using ssh with GitHub and you have gcc compiler in your system. 

1. Clone the source code from this repository

```
git clone https://github.com/Analytics4MD/A4MD_termination_trajectories.git
```

2. Create your conda environment 

```
cd A4MD_termination_trajectories/
conda create --name validation --file environment.yml
conda activate validation
```
The execution of previous commands should create a conda environment that includes all the dependencies that we need to properly run the notebook.

### Testing dependencies
After installing and activating our conda environment, try executing the first two cells in the Jupyter notebook [A4MD_Comparison_early_termination](./A4MD_Comparison_early_termination.ipynb). If no error arises, you are ready to use the Jupyter Notebook. 

## Generating early-terminated trajectories

To generate the early-terminated trajectories you should use our [execute script](./execute.sh). The execute script is in charge of generating the annotations to terminate the trajectories based on some user-defined parameters and then generate the early-terminated trajectories. 

The bash script have some parameters that the user can define:

* `DEBUG` - **"False" or "True" -  Debug mode defines how verbose will the output be**
* `termination_criterion` - **"termination_ess"  or "termination_lev" - Defines what termination criterio will be used in the generation of the early-terminated trajectories**
* `range_min` - **Min LEV range that considers a folded protein**
* `range_max` - **Max LEV range that considers a folded protein**
* `range_tolerance` - **LEV tolerance to remain in range (in %)**
* `window_lev` - **Sliding window size for LEV**
* `stable_th_lev` - **Minimum number of frames considered in window to decide termination for LEV**
* `var_th_ess` - **Minimum average differential variability threshold for ESS (in %)**
* `window_ess` - **Sliding window size for ESS termination**

After defining these values (or using the default), you should be ready to execute this script.

```
./execute.sh
```
The execution of this script will create three important folders:

* `./files/xtc_files`: This is the folder where the full trajectories are downloaded from the [Dataverse](https://doi.org/10.7910/DVN/ML5607)
* `./results/annotations`: This folder has the annotations that are generated to use in the early termination stage
* `./results/trajectories`: This folder has the early-terminated trajectories

After you have the early-terminated trajectories, you can use the [A4MD_Comparison_early_termination](./A4MD_Comparison_early_termination.ipynb) file to visualize the comparison and some other plots (e.g., Free Energy Surface).


## Visual comparison of early terminated trajectories

### Preparing the Notebook

There are many things that can be customized in this Notebook, e.g. the paths to the trajectories (full, trimmed)

* `stride` - **Stride for loading trajectory**
* `selection` - **Used for the RMSD comparison. This option can be "protein" or "all"**
* `input_dirs` - **Path fot the full trajectories of the simulation.**
* `trimmed_trajectories_1` - **Path for the trimmed trajectories of the simulation using LEV**
* `trimmed_trajectories_2` - **Path for the trimmed trajectories of the simulation using ESS**
* `top_file` - **Topology file for the trajectories**
* `nstates` - **Number of states for the MSM model and the PCCA+**
* `anotations`- **Path for the annotations files. The annotations files contains the information of the LEV and ESS termination**
* `SAVE_FRAMES` - **If True, the frames will be saved in the folder frames_closest**
* `frames_closest_folder` - **Folder where the frames will be saved**
* `dist_cmap` - **Color map for the energy plots**
* `size` - **Size of the point in the plots**


## Related Publications

<i class="fa fa-file-text-o"></i> Silvina Caino-Lores, Michel Cuendet, Jack Marquez, Ekaterina Kots, Trilce Estrada, Ewa Deelman, Harel Weinstein, and Michela Taufer.
<b>Runtime steering of molecular dynamics simulations through in situ analysis and annotation of collective variables.</b>
<i>ACM Proceedings of the Platform for Advanced Scientific Computing Conference.</i>
ACM (2023). <a href="https://dl.acm.org/doi/pdf/10.1145/3592979.3593420" target="_blank">[link]</a>

<i class="fa fa-file-text-o"></i> Harshita Sahni, Hector Carrillo-Cabada, Ekaterina Kots, Silvina Caino-Lores, Jack Marquez, Ewa Deelman, Michel Cuendet, Harel Weinstein, Michela Taufer, and Trilce Estrada.
<b>Online Boosted Gaussian Learners for In-Situ Detection and Characterization of Protein Folding States in Molecular Dynamics Simulations.</b>
<i>2023 IEEE 19th International Conference on e-Science (e-Science).</i>
IEEE (2023). <a href="https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=10254895" target="_blank">[link]</a>

<i class="fa fa-file-text-o"></i> Hector Carrillo-Cabada, Jeremy Benson, Asghar Razavi, Brianna Mulligan, Michel A. Cuendet, Harel Weinstein, Michela Taufer, and Trilce Estrada.
<b>A Graphic Encoding Method for Quantitative Classification of Protein Structure and Representation of Conformational Changes</b>
<i>IEEE/ACM Transactions on Computational Biology and Bioinformatics (IEEE/ACM TCBC).</i>
(2020). <a href="https://ieeexplore.ieee.org/document/8859247/" target="_blank">[link]</a>

<i class="fa fa-file-text-o"></i> Tu Mai Anh Do, Loic Pottier, Stephen Thomas, Rafael Ferreira da Silva, Michel A. Cuendet, Harel Weinstein, Trilce Estrada, Michela Taufer, and Ewa Deelman.
<b>A Novel Metric to Evaluate In Situ Workflows</b>
<i>In Proceedings of the International Conference on Computational Science (ICCS), pp. 1 – 14.</i>
(2020). <a href="https://scitech.isi.edu/wordpress/wp-content/papercite-data/pdf/do2020iccs.pdf" target="_blank">[link]</a>

<i class="fa fa-file-text-o"></i> Michela Taufer, Trilce Estrada, and Travis Johnston.
<b>A Survey of Algorithms for Transforming Molecular Dynamics Data into Metadata for In Situ Analytics based on Machine Learning Methods</b>
<i>Issue of Philosophical Transactions A., 378(2166):1-11.</i>
(2020). <a href="https://royalsocietypublishing.org/doi/full/10.1098/rsta.2019.0063" target="_blank">[link]</a>

[More references](https://analytics4md.org/)

## Copyright and License


Copyright (c) 2022, Global Computing Lab

A4MD is distributed under terms of the [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0) with LLVM Exceptions.

See [LICENSE](https://github.com/Analytics4MD/A4MD/blob/PASC/LICENSE) for more details.

## Acknowledgments

This research was supported by the National Science Foundation (NSF) under grant numbers 1741057, 1841758, 2138811, 2223704 and 2331152; 
the Oak Ridge Leadership Computing Facility under allocation CSC427; 
the Extreme Science and Engineering Discovery Environment (XSEDE) under allocation TG-CIS200053;
and IBM through a Shared University Research Award.

## Contact 
Dr. Michela Taufer

[Global Computing Lab](https://globalcomputing.group/)

University of Tennessee

taufer@acm.org

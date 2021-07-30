# Introduction
Computational Reverse Engineering Analysis for Scattering Experiments (CREASE) is a two step approach that combines a genetic algorithm (GA) as the first step with molecular simulations based reconstruction as the second step to determine structural features of assembled structures from aggregate dimensions to the chain- and monomer- level packing for a given input scattering profile. The `crease_ga` package, developed by Michiel Wessels, Zijie Wu, Ziyu Ye, Christian Heil and Arthi Jayaraman at University of Delaware, is an open-source package that implements the GA portion of CREASE.

If using this code, please cite the following references for the method:

1. Beltran-Villegas, D. J.; Wessels, M. G.; Lee, J. Y.; Song, Y.; Wooley, K. L.; Pochan, D. J.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments on Amphiphilic Block Polymer Solutions. J. Am. Chem. Soc. 2019, 141, 14916−14930.

1. Wessels, M. G.; Jayaraman, A. Computational Reverse-Engineering Analysis of Scattering Experiments (CREASE) on Amphiphilic Block Polymer Solutions: Cylindrical and Fibrillar Assembly. Macromolecules 2021, 54, 783-796.

# Installation
To install this package on a linux machine, follow these steps:

1. We will be using Anaconda to handle the python environment. So first, we need to [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) on our machine. We recommend miniconda for faster installation and lower storage space consumption.

1. Download the package. If you have git installed, this can be done using the following command from a terminal:
```
git clone https://github.com/arthijayaraman-lab/crease_ga
```

Or follow the guidance [here](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository-from-github/cloning-a-repository)

2. Create a new conda environment and install all the dependencies.
   - Navigate into the root directory of the cloned package, and create a fresh conda environment with dependencies using
     ```
     conda env create -f environment.yml
     ```
   - Activate the environment with
     ```
     conda activate crease_ga
     ```

3. Install the `crease_ga` package using
   ```
   pip install .
   ```
You are all set to use the `crease_ga` package! Remember to activate the proper environment every time by with `conda activate crease_ga`.

# Getting started
Follow the jupyter notebook tutorials (link) in the `Tutorials` folder, or refer to our documentations at readthedocs.io

# Contact us
If you have any questions or feedbacks, please let us know by emailing to zijiewu@udel.edu

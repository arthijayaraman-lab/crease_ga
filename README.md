#HAHAHA
# Introduction
Computational Reverse Engineering Analysis for Scattering Experiments (CREASE) originally developed by Daniel Beltran-Villegas, Michiel Wessels, and Arthi Jayaraman [Refs 1, 2, and 3 below] is a two step approach that combines a genetic algorithm (GA) as the first step with molecular simulations based reconstruction as the second step to determine structural features of assembled structures from aggregate dimensions to the chain- and monomer- level packing for a given input scattering profile. 

This open source crease_ga package is co-created and maintained by Zijie Wu, Ziyu Ye, Christian Heil and Arthi Jayaraman at University of Delaware. 

If using this code, please cite the following references for the method:

1. Beltran-Villegas, D. J.; Wessels, M. G.; Lee, J. Y.; Song, Y.; Wooley, K. L.; Pochan, D. J.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments on Amphiphilic Block Polymer Solutions. J. Am. Chem. Soc. 2019, 141, 14916−14930. [link to article](https://pubs.acs.org/doi/abs/10.1021/jacs.9b08028)

1. Wessels, M. G.; Jayaraman, A. Computational Reverse-Engineering Analysis of Scattering Experiments (CREASE) on Amphiphilic Block Polymer Solutions: Cylindrical and Fibrillar Assembly. Macromolecules 2021, 54, 783-796. [link to article](https://pubs.acs.org/doi/abs/10.1021/acs.macromol.0c02265)

1. Wessels, M. G.; Jayaraman, A. Machine Learning Enhanced Computational Reverse Engineering Analysis for Scattering Experiments (CREASE) to Determine Structures in Amphiphilic Polymer Solutions. ACS Polym. Au 2021, Advance Article. [link to article](https://pubs.acs.org/doi/abs/10.1021/acspolymersau.1c00015) 

1. Heil, C. M.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments of Assembled Binary Mixture of Nanoparticles. ACS Mater. Au 2021, Advance Article. [link to article](https://pubs.acs.org/doi/10.1021/acsmaterialsau.1c00015) 

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
Follow the jupyter notebook tutorials (link) in the `Tutorials` folder (currently under construction), or refer to our documentations at readthedocs.io (currently under construction).

# Contact us
If you have any questions or feedbacks, please let us know by emailing to zijiewu@udel.edu

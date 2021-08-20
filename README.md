# Introduction
Computational Reverse Engineering Analysis for Scattering Experiments (CREASE) originally developed by Daniel Beltran-Villegas, Michiel Wessels, and Arthi Jayaraman (References 1, 2, and 3 below) is a two step approach that takes as input small angle scattering profile [I(q) vs. q] and provides as output structural information of the assembled structure.  It combines a genetic algorithm (GA) as the first step with molecular simulations based reconstruction as the second step. For macromolecular solutions at dilute concentrations, when the primary component of [I(q) vs. q] is the form factor of the assembled structure, the two steps of CREASE provide structural information about the assembled structure ranging from dimensions of the domains as well as the chain- and monomer- level packing (as shown in References 1, 2, 3, and 5). The work in Reference 4 below describes early implementation of CREASE to analyze input small angle scattering profiles when the structure factor is the primary component of I(q) vs. q. 

This open source crease_ga package has been co-created by Zijie Wu, Ziyu Ye, Christian Heil, and Arthi Jayaraman at University of Delaware.

If you use this code, please cite the one or more of the relevant references from the list below:

1. Original Article on CREASE for spherical micelles:  

   Beltran-Villegas, D. J.; Wessels, M. G.; Lee, J. Y.; Song, Y.; Wooley, K. L.; Pochan, D. J.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments on Amphiphilic Block Polymer Solutions. J. Am. Chem. Soc. 2019, 141, 14916−14930. [link to article](https://pubs.acs.org/doi/abs/10.1021/jacs.9b08028)

2. Extension of CREASE for cylindrical and elliptical micelles: 

   Wessels, M. G.; Jayaraman, A. Computational Reverse-Engineering Analysis of Scattering Experiments (CREASE) on Amphiphilic Block Polymer Solutions: Cylindrical and Fibrillar Assembly. Macromolecules 2021, 54, 783-796. [link to article](https://pubs.acs.org/doi/abs/10.1021/acs.macromol.0c02265)

3. Machine Learning Enhanced CREASE:  

   Wessels, M. G.; Jayaraman, A. Machine Learning Enhanced Computational Reverse Engineering Analysis for Scattering Experiments (CREASE) to Determine Structures in Amphiphilic Polymer Solutions. ACS Polym. Au 2021, Advance Article. [link to article](https://pubs.acs.org/doi/abs/10.1021/acspolymersau.1c00015) 

4. Extension of CREASE's Genetic Algorithm Step to Handle Structure Factors:  

   Heil, C. M.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments of Assembled Binary Mixture of Nanoparticles. ACS Mater. Au 2021, Advance Article. [link to article](https://pubs.acs.org/doi/10.1021/acsmaterialsau.1c00015) 

5. Extension of CREASE for vesicles as well as the the ability to estimate polydispersity in dimensions of the domains in the assembled structure and distribution of molecules between the different domains of the assembled structure: 

   Ye, Z.; Wu, Z.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments(CREASE) on Vesicles Assembled from Amphiphilic Macromolecular Solutions. (under review)

# Installation
To install this package on a linux machine, follow these steps:

1. We will be using Anaconda to handle the python environment. So first, we need to [install anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) on our machine. You can also consider installing [miniconda](https://docs.conda.io/en/latest/miniconda.html) for faster installation and lower storage consumption (Note: Miniconda installation package  does not include the 'Anaconda command prompt'. See step 2 for detail)

1. Afterwards, you will need a bash terminal to work with. If you are on a linux or MacOS machine, you can directly launch the terminal. If you are on a Windows machine, you can consider installing [Windows Subsytem for Linux (WSL)](https://ubuntu.com/wsl), or directly use the 'Anaconda prompt' app that comes with anaconda installation. (should be accessible from the Start Menu if you have completed step 1 correctly).

1. Download the package. If you have git installed, this can be done using the following command from a terminal:
```
git clone https://github.com/arthijayaraman-lab/crease_ga
```

Or download directly from github webpage by following the guidance [here](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository-from-github/cloning-a-repository)

1. Create a new conda environment and install all the dependencies. 
   - Navigate into the root directory of the cloned package, and create a fresh conda environment with dependencies using
     ```
     conda env create -f environment.yml
     ```
   - Activate the environment with
     ```
     conda activate crease_ga
     ```

1. Install the `crease_ga` package using
   ```
   pip install .
   ```
You are all set to use the `crease_ga` package! Remember to activate the proper environment every time by with `conda activate crease_ga`. You can run the Jupyter notebook tutorials by using
    ```
    jupyter notebook
    ```

- **NOTE**: if you intend to run this on a supercomputing cluster, you will need to follow the steps to create a python environment of the corresponding cluster.

# Getting started
Follow the jupyter notebook tutorials (CREASE_vesicles_tutorial-NEW.ipynb) in the `Tutorials` folder, or refer to our documentations at readthedocs.io (currently under construction).

# Contact us
If you have any questions or feedback, please let us know by emailing arthij@udel.edu and zijiewu@udel.edu

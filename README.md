[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/arthijayaraman-lab/crease_ga/master) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Documentation Status](https://readthedocs.org/projects/crease-ga/badge/?version=latest)](https://crease-ga.readthedocs.io/en/latest/?badge=latest) [![DOI](https://zenodo.org/badge/387868834.svg)](https://zenodo.org/badge/latestdoi/387868834)

# Introduction
Computational Reverse Engineering Analysis for Scattering Experiments (CREASE) originally developed by Daniel Beltran-Villegas, Michiel Wessels, and Arthi Jayaraman (References 1, 2, and 3 below) is a two step approach that takes as input small angle scattering profile [I(q) vs. q] and provides as output structural information of the assembled structure.  It combines a genetic algorithm (GA) as the first step with molecular simulations based reconstruction as the second step. For macromolecular solutions at dilute concentrations, when the primary component of [I(q) vs. q] is the form factor of the assembled structure, the two steps of CREASE provide structural information about the assembled structure ranging from dimensions of the domains as well as the chain- and monomer- level packing (as shown in References 1, 2, 3, and 5). The work in Reference 4 below describes early implementation of CREASE to analyze input small angle scattering profiles when the structure factor is the primary component of [I(q) vs. q].

This open source `crease_ga` package has been co-created by Zijie Wu, Ziyu Ye, Christian Heil, and Arthi Jayaraman at University of Delaware.

__If you use this code, please cite one or more of the relevant references from the list below:__

1. Original Article on CREASE for spherical micelles:  

   Beltran-Villegas, D. J.; Wessels, M. G.; Lee, J. Y.; Song, Y.; Wooley, K. L.; Pochan, D. J.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments on Amphiphilic Block Polymer Solutions. J. Am. Chem. Soc. 2019, 141, 14916−14930. [link to article](https://pubs.acs.org/doi/abs/10.1021/jacs.9b08028)

2. Extension of CREASE for cylindrical and elliptical micelles: 

   Wessels, M. G.; Jayaraman, A. Computational Reverse-Engineering Analysis of Scattering Experiments (CREASE) on Amphiphilic Block Polymer Solutions: Cylindrical and Fibrillar Assembly. Macromolecules 2021, 54, 783-796. [link to article](https://pubs.acs.org/doi/abs/10.1021/acs.macromol.0c02265)

3. Machine Learning Enhanced CREASE:  

   Wessels, M. G.; Jayaraman, A. Machine Learning Enhanced Computational Reverse Engineering Analysis for Scattering Experiments (CREASE) to Determine Structures in Amphiphilic Polymer Solutions. ACS Polym. Au 2021, 1, 3, 153-164. [link to article](https://pubs.acs.org/doi/abs/10.1021/acspolymersau.1c00015) 

4. Extension of CREASE's Genetic Algorithm Step to Handle Structure Factors:  

   Heil, C. M.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments of Assembled Binary Mixture of Nanoparticles. ACS Mater. Au 2021, 1, 2, 140-156. [link to article](https://pubs.acs.org/doi/10.1021/acsmaterialsau.1c00015) 

5. Extension of CREASE for vesicles as well as the the ability to estimate polydispersity in dimensions of the domains in the assembled structure and distribution of molecules between the different domains of the assembled structure: 

   Ye, Z.; Wu, Z.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments (CREASE) on Vesicles Assembled from Amphiphilic Macromolecular Solutions. JACS Au 2021, 1, 11, 1925-1936. [link to article](https://pubs.acs.org/doi/10.1021/jacsau.1c00305)

6. Machine Learning Enhanced CREASE for Structure Factors with Expansion to Nanoparticle Solutions:  

   Heil, C. M.; Patil, A.; Dhinojwala, A.; & Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments (CREASE) with Machine Learning Enhancement to Determine Structure of Nanoparticle Mixtures and Solutions. ACS Central Science 2022, 8, 7, 996-1007. [link to article](https://pubs.acs.org/doi/full/10.1021/acscentsci.2c00382) 
   
7. Machine Learning Enhanced CREASE for Semi-flexible Fibrils:  

   Wu, Z. & Jayaraman, A.  Machine learning enhanced computational reverse-engineering analysis for scattering experiments (CREASE) for analyzing fibrillar structures in polymer solutions. Macromolecule 2022, 55, 24, 11076-11091. [link to article](https://pubs-acs-org.udel.idm.oclc.org/doi/full/10.1021/acs.macromol.2c02165)
   
8. Machine Learning Enhanced CREASE for Simultaneous Form Factor and Structure Factor Elucidation for Concentrated Macromolecular Solutions (e.g., micelles, polymer coated nanoparticles):  

   Heil, C. M.; Ma, Y.; Bharti, B.; & Jayaraman, A.  Computational Reverse-Engineering Analysis for Scattering Experiments for Form Factor and Structure Factor Determination ('P(q) and S(q) CREASE'). JACS Au 2023, ASAP. [link to article](https://pubs-acs-org.udel.idm.oclc.org/doi/10.1021/jacsau.2c00697#)

# Installation

To install this package on a Windows, linux, or macOS machine, follow these steps:

1. We will be using Anaconda to handle the python environment. So first, we need to [install anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) on our machine. You can also consider installing [miniconda](https://docs.conda.io/en/latest/miniconda.html) for faster installation and lower storage consumption _only if you are an advanced user_. We recommend installing the full version of anaconda if you have limited experience with python and/or computer programming.

1. At this point, you will need a bash terminal to work with. All the installation steps _after_ this step will need to be accomplished from a terminal, using command line interface. 
    - If you are on a linux or MacOS machine
    
       You can directly launch a terminal.
    - If you are on a Windows machine
    
       If you have installed the full version of Anaconda/Miniconda in Step 1, the most straightforward way to launch a terminal will be using _Anaconda prompt_ that comes with the conda installation. You should be able to find _Anaconda prompt_ from your Start menu. You can also consider installing [Windows Subsytem for Linux (WSL)](https://ubuntu.com/wsl).

1. Download the package. 
   - If you have git installed, this can be done using the following command from a terminal:
     ```
     git clone https://github.com/arthijayaraman-lab/crease_ga
     ```
   - You can also directly download the package as a ZIP file from [our github webpage](https://github.com/arthijayaraman-lab/crease_ga) by following the guidance [here](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository-from-github/cloning-a-repository). If you are following this route, you will need to unzip the package to your desired location after the download.

1. Create a new conda environment and install the package and all its dependencies. 
   - _Navigate into the root directory of the cloned package_. If you are using the anaconda prompt, you can look up [common windows command line prompts](http://www.cs.columbia.edu/~sedwards/classes/2015/1102-fall/Command%20Prompt%20Cheatsheet.pdf). If you are using a unix-based terminal (linux, macOS, or WSL subsystem), you can look up [common commands for Unix](http://www.mathcs.emory.edu/~valerie/courses/fall10/155/resources/unix_cheatsheet.html). Either case, all you will need would probably be _displaying list of files and directories in the current folder_(`dir` for windows, `ls` for unix), and _moving to different directories_(`cd [new_directory_path]` for both windows and unix). You should end up at a directory called `crease_ga` and be able to see files named `setup.py`, `environment.yml` and `README.md` (among others) in the directory.
   - create a fresh conda environment with the package and its dependencies installed using
     ```
     conda env create -f environment.yml
     ```
   - Activate the environment with
     ```
     conda activate crease_ga
     ```

To test if the `crease-ga` package is installed properly, run
```
python3
```
to launch python, and then in the resulting python command line, run
```
import crease_ga
crease_ga.__version__
```

If all installation steps are done properly, you should see the version number of the package printed, and you are all set to use the `crease_ga` package! Remember to activate the proper environment every time by with `conda activate crease_ga`. You can run the Jupyter notebook tutorials with the command
    ```
    jupyter notebook
    ```

- **NOTE1**: If you would like to first try our package by running our tutorial (See #Getting Started), you can directly launch a docker image of our environment to access and interact with our jupyter notebook tutorial from a web browser _without performing any installation steps above_. To do this,  click this badge:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/arthijayaraman-lab/crease_ga/master).

- **NOTE2**: if you intend to run this on a supercomputing cluster, you will need to follow the steps to create a python environment of the corresponding cluster.


# Getting started
Follow the jupyter notebook tutorials (CREASE_vesicles_tutorial-NEW.ipynb) in the `Tutorials` folder, or refer to our documentations at readthedocs.io (currently under construction).

# Contact us
If you have any questions or feedback, please let us know by emailing creasejayaramanlab AT gmail.com.



.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/arthijayaraman-lab/crease_ga/master
   :alt: Binder
 
.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT
 
.. image:: https://zenodo.org/badge/387868834.svg
   :target: https://zenodo.org/badge/latestdoi/387868834
   :alt: DOI


Introduction
============

Computational Reverse Engineering Analysis for Scattering Experiments (CREASE) originally developed by Daniel Beltran-Villegas, Michiel Wessels, and Arthi Jayaraman (References 1, 2, and 3 below) is a two step approach that takes as input small angle scattering profile [I(q) vs. q] and provides as output structural information of the assembled structure.  It combines a genetic algorithm (GA) as the first step with molecular simulations based reconstruction as the second step. For macromolecular solutions at dilute concentrations, when the primary component of [I(q) vs. q] is the form factor of the assembled structure, the two steps of CREASE provide structural information about the assembled structure ranging from dimensions of the domains as well as the chain- and monomer- level packing (as shown in References 1, 2, 3, and 5). The work in Reference 4 below describes early implementation of CREASE to analyze input small angle scattering profiles when the structure factor is the primary component of I(q) vs. q. 

This open source crease_ga package has been co-created by Zijie Wu, Ziyu Ye, Christian Heil, and Arthi Jayaraman at University of Delaware.

References
_________

**If you use this code, please cite one or more of the relevant references from the list below:**


#. 
   Original Article on CREASE for spherical micelles:  

   Beltran-Villegas, D. J.; Wessels, M. G.; Lee, J. Y.; Song, Y.; Wooley, K. L.; Pochan, D. J.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments on Amphiphilic Block Polymer Solutions. J. Am. Chem. Soc. 2019, 141, 14916−14930. `link to article <https://pubs.acs.org/doi/abs/10.1021/jacs.9b08028>`_

#. 
   Extension of CREASE for cylindrical and elliptical micelles: 

   Wessels, M. G.; Jayaraman, A. Computational Reverse-Engineering Analysis of Scattering Experiments (CREASE) on Amphiphilic Block Polymer Solutions: Cylindrical and Fibrillar Assembly. Macromolecules 2021, 54, 783-796. `link to article <https://pubs.acs.org/doi/abs/10.1021/acs.macromol.0c02265>`_

#. 
   Machine Learning Enhanced CREASE:  

   Wessels, M. G.; Jayaraman, A. Machine Learning Enhanced Computational Reverse Engineering Analysis for Scattering Experiments (CREASE) to Determine Structures in Amphiphilic Polymer Solutions. ACS Polym. Au 2021, Advance Article. `link to article <https://pubs.acs.org/doi/abs/10.1021/acspolymersau.1c00015>`_ 

#. 
   Extension of CREASE's Genetic Algorithm Step to Handle Structure Factors:  

   Heil, C. M.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments of Assembled Binary Mixture of Nanoparticles. ACS Mater. Au 2021, Advance Article. `link to article <https://pubs.acs.org/doi/10.1021/acsmaterialsau.1c00015>`_ 

#. 
   Extension of CREASE for vesicles as well as the the ability to estimate polydispersity in dimensions of the domains in the assembled structure and distribution of molecules between the different domains of the assembled structure: 

   Ye, Z.; Wu, Z.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments(CREASE) on Vesicles Assembled from Amphiphilic Macromolecular Solutions. (under review)

.. toctree::
    :caption: Getting started
    :maxdepth: 2

    installation
    tutorials
.. toctree::    
    :caption: Topics
    :maxdepth: 2
    
    topics/shape
.. toctree::    
    :caption: Documentations
    :maxdepth: 1

    documentations/Model

Contact us
__________

If you have any questions or feedback, please let us know by emailing arthij@udel.edu and zijiewu@udel.edu
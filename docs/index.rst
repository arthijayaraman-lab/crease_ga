

.. image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/arthijayaraman-lab/crease_ga/master
   :alt: Binder
 
.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License: MIT
 
.. image:: https://zenodo.org/badge/387868834.svg
   :target: https://zenodo.org/badge/latestdoi/387868834
   :alt: DOI


What is CREASE?
===============

Computational Reverse Engineering Analysis for Scattering Experiments (CREASE) originally developed by Daniel Beltran-Villegas, Michiel Wessels, and Arthi Jayaraman (References 1, 2, and 3 below) is a two step approach that takes as input small angle scattering profile [I(q) vs. q] and provides as output structural information of the assembled structure.  It combines a genetic algorithm (GA) as the first step with molecular simulations based reconstruction as the second step. For macromolecular solutions at dilute concentrations, when the primary component of [I(q) vs. q] is the form factor of the assembled structure, the two steps of CREASE provide structural information about the assembled structure ranging from dimensions of the domains as well as the chain- and monomer- level packing (as shown in References 1, 2, 3, and 5). The work in Reference 4 below describes early implementation of CREASE to analyze input small angle scattering profiles when the structure factor is the primary component of [I(q) vs. q]. 

This open source crease_ga package has been co-created by Zijie Wu, Ziyu Ye, Christian Heil, and Arthi Jayaraman at the University of Delaware and focuses on the first genetic algorithm (GA) step of CREASE.

References
__________

**If you use this code, please cite one or more of the relevant references from the list below:**


#. 
   Original Article on CREASE for spherical micelles:  

   Beltran-Villegas, D. J.; Wessels, M. G.; Lee, J. Y.; Song, Y.; Wooley, K. L.; Pochan, D. J.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments on Amphiphilic Block Polymer Solutions. J. Am. Chem. Soc. 2019, 141, 14916−14930. `link to article <https://pubs.acs.org/doi/abs/10.1021/jacs.9b08028>`_


#. 
   Extension of CREASE for cylindrical and elliptical micelles: 

   Wessels, M. G.; Jayaraman, A. Computational Reverse-Engineering Analysis of Scattering Experiments (CREASE) on Amphiphilic Block Polymer Solutions: Cylindrical and Fibrillar Assembly. Macromolecules 2021, 54, 783-796. `link to article <https://pubs.acs.org/doi/abs/10.1021/acs.macromol.0c02265>`_


#. 
   Machine Learning Enhanced CREASE:  

   Wessels, M. G.; Jayaraman, A. Machine Learning Enhanced Computational Reverse Engineering Analysis for Scattering Experiments (CREASE) to Determine Structures in Amphiphilic Polymer Solutions. ACS Polym. Au 2021, 1, 3, 153-164. `link to article <https://pubs.acs.org/doi/abs/10.1021/acspolymersau.1c00015>`_ 


#. 
   Extension of CREASE's Genetic Algorithm Step to Handle Structure Factors:  

   Heil, C. M.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments of Assembled Binary Mixture of Nanoparticles. ACS Mater. Au 2021, 1, 2, 140-156. `link to article <https://pubs.acs.org/doi/10.1021/acsmaterialsau.1c00015>`_ 


#. 
   Extension of CREASE for vesicles as well as the the ability to estimate polydispersity in dimensions of the domains in the assembled structure and distribution of molecules between the different domains of the assembled structure: 

   Ye, Z.; Wu, Z.; Jayaraman, A. Computational Reverse-Engineering Analysis for Scattering Experiments (CREASE) on Vesicles Assembled from Amphiphilic Macromolecular Solutions. JACS Au 2021, 1, 11, 1925-1936. `link to article <https://pubs.acs.org/doi/10.1021/jacsau.1c00305>`_


#.
   Machine Learning Enhanced CREASE for Structure Factors with Expansion to Nanoparticle Solutions:

   Heil, C. M.; Patil, A.; Dhinojwala, A.; & Jayaraman, A. (2022). Computational Reverse-Engineering Analysis for Scattering Experiments (CREASE) with Machine Learning Enhancement to Determine Structure of Nanoparticle Mixtures and Solutions. ACS Central Science. `link to article <https://pubs-acs-org.udel.idm.oclc.org/doi/10.1021/acscentsci.2c00382>`_


#. 
   Machine Learning Enhanced CREASE for Semi-flexible Fibrils:  

   Wu, Z. & Jayaraman, A. (2022). Machine learning enhanced computational reverse-engineering analysis for scattering experiments (CREASE) for analyzing fibrillar structures in polymer solutions. Macromolecule. `link to article <https://pubs-acs-org.udel.idm.oclc.org/doi/full/10.1021/acs.macromol.2c02165>`_ 
   
#.
   Machine Learning Enhanced CREASE for Simultaneous Form Factor and Structure Factor Elucidation for Concentrated Macromolecular Solutions (e.g., micelles, polymer coated nanoparticles):  

   Heil, C. M.; Ma, Y.; Bharti, B.; & Jayaraman, A. (ArXiv preprint). Computational Reverse-Engineering Analysis for Scattering Experiments for Form Factor and Structure Factor Determination ('P(q) and S(q) CREASE'). `link to article <https://arxiv.org/abs/2212.03154>`_ 


Contact us
__________

If you have any questions or feedback, please let us know by emailing creasejayaramanlab AT gmail.com.

.. toctree::
    :caption: Getting started
    :maxdepth: 1
    
    self
    installation

.. toctree::
    :caption: Tutorials
    :maxdepth: 1

    CREASE_GA Vesicle <tutorials>

.. toctree::
    :caption: Crease_ga Workflow
    :maxdepth: 1

    example_workflow

.. toctree::    
    :caption: Manual for Codes
    :maxdepth: 1

    documentations/Model
    documentations/shape-vesicle
    documentations/shape-micelle
    documentations/adaptation_params
    documentations/utils

.. toctree::    
    :caption: Miscellaneous
    :maxdepth: 1
    
    Frequently Asked Questions <topics/faq>
    How can you contribute? <topics/plugin>
    topics/feedback


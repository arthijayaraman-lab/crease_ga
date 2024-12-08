Case Study I: Analyzing Small Angle Scattering Profiles (1D) from Soft Materials with Dispersity in Particle Size and Shapes using CREASE
=======================================================================================================================================

In this case study, we use the CREASE method to analyze the dispersity in both the size and shape of nanoparticles from their 1D scattering profile. A computational method that seeks to analyze the dispersity in both the size and shape simultaneously must have the ability to output multiple solutions. This has traditionally been a challenge for analytical models. The CREASE genetic algorithm outputs multiple solutions to an input scattering profile and ranks them according to their fitness. Hence it is capable of analyzing the dispersity in both the size and shape of nanoparticles simultaneously. Below we walk you through the 4 steps involved in training a Machine Learning (ML) model that links the structural features of a nanoparticle system (with dispersity in both the size and shape) directly to its 1D scattering profile. Once the ML model is trained, it can be incorporated into the CREASE loop. In this case, CREASE takes a 1D scattering profile as an input and outputs multiple sets of structural features whose computed scattering profiles closely match the input scattering profile.        

Step 1: Identifying the Structural Features of the System  
----------------------------------------------------------

The first step in training an ML model (that is to be incorporated in the CREASE loop) is the identification of structural features that are relevant to the system. In this study, the size of a nanoparticle is defined by its volumetric radius. Therefore, the structural features that capture the size distribution of the nanoparticle system are the mean and standard deviation of the volumetric radii. The shape of a nanoparticle is defined by its aspect ratio in this study. Therefore, the structural features that capture the shape distribution of the nanoparticle system are the mean and standard deviation of the aspect ratios. We define an additional structural feature, volume fraction, that defines the degree of packing of nanoparticles. **Figure 1** provides a pictorial description of the structural features.   

.. figure:: CasestudyI_Step_1.png
   :class: with-border

   Figure 1.: Definitions of structural features relevant to this case study. Structural features are a set of physical descriptors that fully define a 3 dimensional (3D) structure of the system. For a system of nanoparticles with dispersity in both size and shape, we identify a set of five structural features as shown above.    

During this step it is important to spend some time and think about the expected output from CREASE. The final output from CREASE is going to be multiple sets of structural features that are identified here.         

Step 2:	Generating 3D Structures for Varying Values of Structural Features
----------------------------------------------------------------------------

.. figure:: CasestudyI_Step_2.png
   :class: with-border

   Figure 2.: Shows four 3 dimensional (3D) structures generated using the CASGAP **[1]** program. The 4 sets of structural features that were input to CASGAP to obtain the 3D structures are also shown. 


3.	Calculating Scattering Profiles for the 3D Structures Generated
---------------------------------------------------------------------



.. figure:: CasestudyI_Step_3.png
   :class: with-border 

Figure 3.: .

4.	How has CREASE been used so far?
----------------------------------------

CREASE method has been used to interpret small angle scattering results to 

a. *Identify relevant dimensions of assembled structures in polymer solutions at dilute concentrations* **[5-9]**: CREASE has  been applied to characterize structure of the ‘primary particle’ using scattering profiles I(q) ~ P(q) (*i.e.*, conditions where S(q) is ~1) for a variety of ‘primary particles’ (micelles **[6, 7, 9]**, vesicles **[8]**, and fibrils **[5]**) bypassing the need for an analytical model. 

b.	*Understand the amorphous structure of spherical particles at high concentrations regardless of extent of mixing/segregation*: CREASE has also been extended to analyze S(q) part of the scattering profiles from concentrated binary mixture of polydisperse spherical nanoparticles (i.e., P(q) is a sphere form factor) to determine the extent of segregation/mixing of the two types of nanoparticles and the precise mixture composition **[4, 10]**. 

c.	*Elucidate the amorphous structure of particles / micelles in solutions, with unknown primary particle form and unknown assembled/dispersed structure* **[11]**: Most recently, for systems where one does not know the P(q) or S(q) a priori, CREASE has been extended to simultaneously interpret structural information held in P(q) and S(q) and appropriately called ‘P(q) and S(q) CREASE’ **[11]**.

*CREASE has taken as input 1D SAXS profiles and/or SANS profiles*: In the studies above, the input to CREASE has been (i) a single SAXS profile of the system, or (ii) one SAXS profile and a one SANS profile of the same system, or (iii) multiple SANS profiles with contrast matching one or the other component(s) in the system with the solvent. Next development steps of CREASE development are focused on 2D profiles for soft materials that show anisotropy in the assembled structure.

*CREASE with Debye method vs. ML-model for computed scattering profile calculation*: In earlier implementations of CREASE, the Debye method for computed scattering profile calculation was used; as noted above this calculation was initially found to be quite time consuming. In following work, the structure generation (done in every step of Debye method) was found to more computationally intensive while the computed scattering calculations using Debye method have been made faster than in previous implementations. The machine learning (ML) enhanced CREASE-GA, with a well-trained ML model avoids both Debye equation based computed scattering calculation and the three-dimensional real space structure generation in the optimization loop, making is significantly faster than using Debye method (*e.g.*, one can complete CREASE-GA optimization is less than an hour on a laptop with a pre-trained ML model!)


References
__________

#.
   Brisard, S.; Levitz, P., *Small-angle scattering of dense, polydisperse granular porous media: Computation free of size effects.*
   **Phys. Rev. E 2013, 87 (1), 013305.** (`link <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.87.013305>`_)

#.
   Olds, D. P.; Duxbury, P. M., *Efficient algorithms for calculating small-angle scattering from large model structures.*
   **Journal of Applied Crystallography 2014, 47 (3), 1077-1086.** (`link <https://journals.iucr.org/j/issues/2014/03/00/kk5148/index.html>`_)

#.
   Schmidt-Rohr, K., *Simulation of small-angle scattering curves by numerical Fourier transformation.*
   **Journal of Applied Crystallography 2007, 40 (1), 16-25.** (`link <https://onlinelibrary.wiley.com/iucr/doi/10.1107/S002188980604550X>`_)

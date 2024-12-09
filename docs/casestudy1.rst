Case Study I: Analyzing Small Angle Scattering Profiles (1D) from Soft Materials with Dispersity in Particle Size and Shapes using CREASE
=======================================================================================================================================

In this case study, we use the CREASE method to analyze the dispersity in both the size and shape of nanoparticles from their 1D scattering profile. A computational method that seeks to analyze the dispersity in both the size and shape simultaneously must have the ability to output multiple solutions. This has traditionally been a challenge for analytical models. The CREASE genetic algorithm outputs multiple solutions to an input scattering profile and ranks them according to their fitness. Hence it is capable of analyzing the dispersity in both the size and shape of nanoparticles simultaneously. Below we walk you through the 4 steps involved in training a Machine Learning (ML) model that links the structural features of a nanoparticle system (with dispersity in both the size and shape) directly to its 1D scattering profile. Once the ML model is trained, it can be incorporated into the CREASE loop. In this case, CREASE takes a 1D scattering profile as an input and outputs multiple sets of structural features whose computed scattering profiles closely match the input scattering profile.        

Step 1: Identifying the Structural Features of the System  
----------------------------------------------------------

The first step in training an ML model (that is to be incorporated in the CREASE loop) is the identification of structural features that are relevant to the system. In this study, the size of a nanoparticle is defined by its volumetric radius. Therefore, the structural features that capture the size distribution of the nanoparticle system are the mean and standard deviation of the volumetric radii. The shape of a nanoparticle is defined by its aspect ratio in this study. Therefore, the structural features that capture the shape distribution of the nanoparticle system are the mean and standard deviation of the aspect ratios. We define an additional structural feature, volume fraction, that defines the degree of packing of nanoparticles. **Figure 1** provides a pictorial description of the structural features identified for this study.   

.. figure:: CasestudyI_Step_1.png
   :class: with-border

   Figure 1.: Definitions of structural features relevant to this case study. Structural features are a set of physical descriptors that fully define a 3 dimensional (3D) structure of the system. For a system of nanoparticles with dispersity in both size and shape, we identify a set of five structural features as shown above.    

During this step it is important to spend some time and think about the expected output from CREASE. The final output from CREASE is going to be multiple sets of structural features that are identified here.         

Step 2:	Generating 3D Structures for Varying Values of Structural Features
----------------------------------------------------------------------------
In step 2, we generate 3000 3D structures by varying the structural features defined in the previous step. For this study, the mean volumetric radius of the system was varied uniformly in the range of 20-100 Angstroms. The standard deviation of the volumetric radii was varied uniformly to lie between 0-50% of the mean volumetric radii. The mean aspect ratio of the system was varied uniformly between 0.8-4.0. The standard deviation of the aspect ratio was between 0-50% of the mean aspect ratio. Volume fraction was uniformly varied in the range of [0.05, 0.15]. The range of variation of structural features is determined by manual matching. Manual matching is (also known as sensitivity analysis) carried out before a large scale generation of 3D structures to ensure that the features present in experimental scattering profiles are consistent with the features in the computed scattering profile. Manual matching is also necessary to ensure that the structural features in the chosen range have a noticeable effect on the computed scattering profiles. 4 out of the 3000 sets of structural features and their corresponding 3D representations are shown in **Figure 2**.   

.. figure:: CasestudyI_Step_2.png
   :class: with-border

   Figure 2.: Shows four 3 dimensional (3D) structures generated using the CASGAP **[1]** program. The 4 sets of structural features that were input to CASGAP to obtain the 3D structures are also shown. 

To generate representative 3D structures by using sets of structural features as input, we use the CASGAP **[1]** program developed in the Jayaraman lab. CASGAP generates 3D structures of nanoparticles for a user defined distribution of particle size and shapes. The orientation of the nanoparticles can also by controlled by the user with a kappa (anisotropy) parameter. In this study the kappa parameter was set to 0 for all 3D structures which indicates completely disordered nanoparticle orientations.  

Step 3:	Calculating Scattering Profiles for the 3D Structures Generated
---------------------------------------------------------------------



.. figure:: CasestudyI_Step_3.png
   :class: with-border 

Figure 3.: Shows four computed scattering profiles and the corresponding 3 dimensional (3D) structures. The scattering profiles were computed from the 3D structures using a physics based equation. 

Step 4.	Training a Machine Learning Model that Directly Links Structural Features to the Computed Scattering Profiles
----------------------------------------



.. figure:: CasestudyI_Step_4.png
   :class: with-border 

Figure 4.: Graphical representation of training an XGBoost Machine Learning (ML) model to directly link the structural features of a nanoparticle system to its computed scattering profile. 80% of the scattering profiles computed in step 3 are selected randomly and used as a training dataset for the ML model. The predictions of the ML model are validated by using the remaining 20% of the dataset (test dataset). 

Incorporating the Trained ML Model in CREASE to Analyze the Dispersity in the Size and Shapes of Nanoparticles from their Experimental Scattering Profile
----------------------------------------

References
__________

#.
   Brisard, S.; Levitz, P., *Small-angle scattering of dense, polydisperse granular porous media: Computation free of size effects.*
   **Phys. Rev. E 2013, 87 (1), 013305.** (`link <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.87.013305>`_)


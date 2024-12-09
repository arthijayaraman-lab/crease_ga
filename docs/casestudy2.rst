Case Study II: Analysis of in silico 2D scattering profiles with distributions of particle shapes, sizes and orientational order using CREASE-2D
=======================================================================================================================================

        

Step 1: Structural Feature Identification And Structure Generation
------------------------------------------------------------------

The in-silico system of spheroidal particles has been characterized by 6 structural features, namely :math:`R_\mu`, R<sub>σ</sub>, γ<sub>μ</sub>, γ<sub>σ</sub>, κ and φ. We used the CASGAP method to generate these structures, as it was designed for this type of application. Some relevant details about the identified structural features are provided below and more detailed information can be found in references [1] and [2].

1. R<sub>μ</sub> and R<sub>σ</sub> represent the mean and standard deviation in the volumeteric radius of the spheroidal particles (a measure of their size) following a log-normal distribution.

2. γ<sub>μ</sub> and γ<sub>σ</sub> represent the mean and standard deviation in the spheroidal aspect ratio of the particles, a parameter that directly influences the shape of the particles, also following a log-normal distribution.

3. The orientation of particles is specified by the 3D vector pointing along the major axis of the spheroid. The distribution of orientation follows the 3D von Mises–Fisher (vMF) distribution. Here the mean orientation has been fixed to one of the axes in the laboratory frame (as explained in the manuscript) and the extent of the orientational order is characterized by the κ parameter, which is a measure of the inverse dispersity in orientation. κ = 0 indicates complete lack of the orientational order (isotropic) and κ → ∞ indicates the perfect orientational order (highly anisotropic). 

4.	The concentration of particles is quantified by the volume fraction of particles, ϕ. This parameter is not independently controlled and is evaluated after the CASGAP method generates a structure with the specified distributions of shape, size and orientational order.


.. figure:: case_study_2_files/Figure2_Step1.png
   :class: with-border

   Figure 2.:       

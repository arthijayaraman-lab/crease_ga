Case Study II: Analysis of 2D scattering profiles for *in silico* system of spheroidal particles with dispersities in shapes, sizes and orientational order using CREASE-2D
===========================================================================================================================================================================
This case study demonstrates the implementation of CREASE-2D workflow for an *in silico* system of spheroidal particles as also discussed in the original manuscript for CREASE-2D [1]. The CREASE-2D workflow


Step 1: Structural Feature Identification And Structure Generation
------------------------------------------------------------------

The *in silico* system of spheroidal particles has been characterized by 6 structural features, namely :math:`R_\mu`, :math:`R_\sigma`, :math:`\gamma_\mu`, :math:`\gamma_\sigma`, :math:`\kappa` and :math:`\phi`. Some relevant details about the identified structural features are provided below and more detailed information can be found in references [1] and [2].

#. :math:`R_\mu` and :math:`R_\sigma` represent the mean and standard deviation in the volumeteric radius of the spheroidal particles (a measure of their size) following a log-normal distribution.

#. :math:`\gamma_\mu` and :math:`\gamma_\sigma` represent the mean and standard deviation in the spheroidal aspect ratio of the particles, a parameter that directly influences the shape of the particles, also following a log-normal distribution.

#. The orientation of particles is specified by the 3D vector pointing along the major axis of the spheroid. The distribution of orientation follows the 3D von Mises–Fisher (vMF) distribution. Here the mean orientation has been fixed to one of the axes in the laboratory frame (as explained in the manuscript) and the extent of the orientational order is characterized by the :math:`\kappa` parameter, which is a measure of the inverse dispersity in orientation. :math:`\kappa=0` indicates complete lack of the orientational order (isotropic) and :math:`\kappa\rightarrow\infty` indicates perfect orientational order (highly anisotropic). 

#. The concentration of particles is quantified by the volume fraction of particles, :math:`\phi`. This parameter is not independently controlled and is evaluated after the CASGAP method generates a structure with the specified distributions of shape, size and orientational order.

Using the defined structural features a dataset of 3000 structural features was first generated. We used the CASGAP method [2] to generate these structures, as it was developed primarily to generate overlap-free configurations of spheroidal particles with pre-defined distributions of shapes, sizes and orientational order. The identified structural features and some representative structures are depicted in Figure 2.

.. figure:: case_study_2_files/Figure2_Step1.png
   :class: with-border
   :width: 900px
   :align: center

   Figure 2.: **(A)** Identified structural features for the *in silico* system. **(B-D)** Representative snapshots of 3D structures displaying variations in size, shape and orientational order, respectively. Figure adapted from reference [1].

**Important Note: CREASE-2D has now been extended to work with structures that are entirely defined by point scatterers, which makes it adaptable to any structural configuration (not just spheroids). Thus Step 1 can be adapted to include any system with defined or identified structural features, and any computational method that generates a 3D structure (filled with point scatterers) that can be manipulated by those defined or identified structural features.** 


Step 2:	Calculation of Scattering Profiles
------------------------------------------

For each of the generated structures in Step 1, 2D scattering intensity :math:`I(q,\theta)` is computed by first computing the scattering amplitude :math:`A(q,\theta)`. Calculation of scattering amplitude involves a single summation term over the entire list of scatterers, which is easier to parallelize over multiple cpus or gpus. The results of such calculations for the current *in silico* system are shown in Figure 3.

.. figure:: case_study_2_files/Figure3_Step2.png
   :class: with-border
   :width: 900px
   :align: center

   Figure 3.: Calculated 2D scattering profiles in cartesian (center) and polar (right) form for representative structures (left) shown for a few samples. Figure adapted from reference [1].

The cartesian form of the 2D scattering intensity :math:`I(q,\theta)` is directly used along with the structural features to obtain the dataset for ML training in Step 3.

Step 3:	Training of Surrogate Machine Learning Model to Predict Scattering Profiles from Structural Features
------------------------------------------------------------------------------------------------------------

The data set of 3000 2D scattering profiles and their corresponding structural features is first split such that 80% of the data (2400 structures) is used for training the ML model and the remaining 20% (600 structures) is used for testing/validation of the ML model’s performance. Currently CREASE-2D implementation uses XGBoost as the ML model due to its exceptional performance and lower scope of overfitting. To use XGBoost, the training data set is reformatted into a table, with each row containing all 6 structural features as well as, three new fields corresponding to :math:`q`, :math:`\theta` and :math:`I(q,\theta)`. The last three fields can be obtained by serializing the cartesian form of the 2D scattering profiles, after appropriate subsampling (to avoid excessive data for efficient memory usage; please see main manuscript [1] for more details).

.. figure:: case_study_2_files/Figure4_Step3.png
   :class: with-border
   :width: 900px
   :align: center

   Figure 4.: **(A)** Learning curve during training of XGBoost model, using R\ :sup:`2` error for both the training (black) and validation (green) data entries. **(B)** Performance of the XGBoost model using the R\ :sup:`2` and the structural similarity index measure (SSIM) scores for all 3000 samples in the data set. (C,D) Original and predicted scattering profiles for a selected few samples from the validation data set, each marked with their R\ :sup:`2` and SSIM scores. Figure adapted from reference [1].

Step 4:	Incorporating the Surrogate ML Model within the Genetic Algorithm (GA) Optimization Loop to Complete CREASE-2D Workflow
-------------------------------------------------------------------------------------------------------------------------------
The final step in the CREASE-2D implementation is to put together the predictive capacity and the speed of the surrogate ML model within the GA optimization loop. Consequently, the 6 structural features are represented as 6 corresponding "genes", which are additionally normalized to the interval 0-1. The input to the GA is an *in silico* "experimental" 2D scattering profile (:math:`I_{exp}(q,\theta)`), which is compared to the ML predicted or "computed" 2D scattering profile (:math:`I_{comp}(q,\theta)`).

.. figure:: case_study_2_files/Figure5_Step4.png
   :class: with-border
   :width: 900px
   :align: center

   Figure 5.: **(A)** Identified structural features for the *in silico* system. **(B-D)** Representative snapshots of 3D structures displaying variations in size, shape and orientational order, respectively. Figure adapted from reference [1].

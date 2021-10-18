Frequently Asked Questions
==========================

.. _section-shape:

How does crease_ga support different assembled shapes and morphologies (e.g., vesicle, micelles, fibrils, etc.)?
________________________________________________________________________________________________________________
And
___

Do I have to use the Debye scattering equation to evaluate a scattering profile (I(q))? Can I use known analytical models, trained neural network models, etc.?
_________________________________________________________________________________________________________________________________________________________

With the crease_ga package, we aim to provide a general framework for analyzing scattering profiles using a genetic algorithm (GA) as the optimization method.

#.
        The GA method can be used in the manner described in the Tutorial to analyze and interpret an input scattering profile by assuming a morphology, defining parameters relevant for that morphology (e.g., dimensions, extent of segregation of molecules in various domains within the morphology), placing scatterers within that assumed morphology, calculating the computational scattering profile using the scatterer placements, and finding the values of the parameters whose computed scattering profile closely matches the input scattering profile.
#.
        The GA method can also be used to fit an analytical model for that assumed morphology; in this case it would be optimizing the parameters defined in the analytical model.
        
In the above, the first approach relies on defining a certain genome (a set of parameters) for the assumed morphology of interest (vesicles, micelles, fibrils, etc.) The detailed implementation of this step and its variations (e.g., incorporation of dispersity, parallelization, etc.) all fall under “shapes”.

Currently, we support two built-in “shapes” in the crease_ga codebase: the vesicle model as implemented in `Ye, Wu and Jayaraman <https://pubs.acs.org/doi/10.1021/jacsau.1c00305>`_, and the micelle model as implemented in `Beltran-Villegas et al.  <https://pubs.acs.org/doi/abs/10.1021/jacs.9b08028>`_
. Both “shapes” are copied as-is from the articles, without any parallelization, to facilitate the usage of these “shapes” in a jupyter notebook for pedagogical purposes. While we plan to expand the codebase to include other “shapes” we have developed and are developing, we also encourage users to implement their own shapes of interest.

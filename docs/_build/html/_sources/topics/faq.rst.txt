Frequently Asked Questions
==========================

.. _section-shape:

How does crease_ga support different chemistries (vesicle, micelles, fibrils, etc.)?
___________________________________________________________________________________

Do I have to use Debye scattering equation to evaluate scattering profile (I(q))? Can I use known analytical models, trained neural network models, etc.?
_________________________________________________________________________________________________________________________________________________________

With the crease_ga package, we aim to provide a general framework for analyzing scattering profiles using genetic algorithm (GA). The crease_ga-vesicle tutorial provides a example implementation of this workflow, which consists of two parts: 

#.
        A method to evaluate scattering profile from a set of input parameters is defined (in the vesicle example, placing scatterers based on defined dimensions, and calculating scattering profile using Debye scattering equation)

#.
        Genetic algorithm (GA) for iteratively evaluating scattering profile from a variety of candidates using algorithm defined in part 1, and searching the "best" set of input parameters that produce the most similar scattering profile with the user-provided target.

The second part is largely universal for all possible crease_ga calculations, and the first part, which can be broadly defined as a process of “obtaining a scattering profile from a certain genome (a set of parameters)", differs by the chemistry of interest (vesicelles, micelles, fibrils, etc.), the chosen algorithm to evaluate I(q) (SASView-like analytical models, or debye scattering equation as in the vesicles example, or a pretrained neural network, etc.), and other details (parallelization, assumed distribution of dispersity, etc.). The detailed implementation of this step is delegated to a collection of implementations called "shapes".

**NOTE**: In this context, "shapes" not only handle difference in the assumed morphology of the chemistry, but also all the other aformentioned details, such as algorithms, parallelization implementation, and so on. So for example, a method to calculate I(q) of vesicles using Debye equation, and a method also to calculate I(q) of vesicles, but using some analytical models, should be considered as two different "shapes". In a more extreme case, two methods to calculate I(q) of fibrils using Debye equation, but one assuming dispersity in fibril radius following Gaussian distributon, and the other assuming dispersity in fibril radius following Schultz distribution, should alos be ocnsdered as two different "shapes", as these different assumptions lead to different implementation along the pathway from input parameters to computed scattering profile.

Each shape requires a set of "shape-specific descriptors", and results in GA predicting a set of "input parameters". "Shape-specific descriptors" are values needed to specifically define some additional dimensions or properties of the shape. These descriptors are defined before a crease-ga run, **cannot be changed during the run and are not predicted by GA". In the case of vesicle, these descriptors include number of scatterers used to describe a chain; the ratio of number of monomers of chemistry A: chemistry B; diameter of monomer of chemistry A and B, *etc.* "Input parameters" are the parameters a user would like to fit a scattering profile to - those that will be solved and predicted by crease-ga. In the case of vesicle, the input parameters include the radius of the vesicle core, thickness of the inner, middle and outer layer, dispersity of core radius, *etc.*.

The full list of shape-specific descriptors and input parameters for vesicles, along with the default values of the former and the default minimum/maximum boundaries for the latter when not specified by the user, can be found at :doc:`the API for shape vesicle. <../documentations/shape-vesicle>`

Currently, we only support two builtin "shapes" in crease_ga codebase: the vesicle model as implemented in Ye, Wu and Jayaraman (under review), and the micelle model as implemented in `Beltran-Villegas et al.  <https://pubs.acs.org/doi/abs/10.1021/jacs.9b08028>`_
. Both "shapes" are copied as-is from the articles, but without any parallelization, to facilitate the usage of these "shapes" in jupyter notebook for pedagogical purposes. While we plan to expand the codebase to include other "shapes" we have developed and are developing, we also encourage users to implement their own shapes of interest. 

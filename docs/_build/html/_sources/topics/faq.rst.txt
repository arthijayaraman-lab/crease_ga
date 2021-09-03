Frequently Asked Questions
==========================

How does crease_ga support different chemistries (vesicle, micelles, fibrils, etc.)?
___________________________________________________________________________________

Do I have to use Debye scattering equation to evaluate scattering profile (I(q))? Can I use known analytical models, trained neural network models, etc.?
_________________________________________________________________________________________________________________________________________________________

With the crease_ga package, we aim to provide a general framework for analyzing scattering profiles using genetic algorithm (GA). The crease_ga-vesicle tutorial provides a example implementation of this workflow, which consists of two parts: 

#.
        A method to evaluate scattering profile from a set of input parameters is defined (in the vesicle example, placing scatterers based on defined dimensions, and calculating scattering profile using Debye scattering equation)

#.
        Genetic algorithm (GA) for iteratively searching the "best" set of input parameters that produce the most similar scattering profile with the user-provided target.


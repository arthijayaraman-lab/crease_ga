Simplified Workflow
===================

Crease_ga allows user to fit scattering profiles within several lines of python codes in an intuitive way. Below is a simplified example of how a user can initiate a new fitting model, read the scattering profile from a file, load a :ref:`section-shape` that describes the expected overall morphology (vesicle, micelles, etc.), the algorithm to calculate scattering profile from a set of dimensions (through placement of scatteres and debye scattering equation, or analytical models, or pre-trained neural network, etc.) and other details, and solve the scattering profile based on the chosen shape.

.. code:: ipython3

   import crease_ga as cga
   #Initialize a model
   m = cga.Model(pop_number = 5, generations = 5, nloci = 7)
   #load a shape    
   m.load_shape(shape='vesicle',shape_params=[24,54,0.5,50.4,50.4,0.55,7],
                                     minvalu = (50, 30, 30, 30, 0.1, 0.0, 0.1),
                                     maxvalu = (400, 200, 200, 200, 0.45, 0.45, 4))
   #Read Iexp(q) from a file                                  
   m.load_iq('../IEXP_DATA/Itot_disper_10_Ain12_B6_Aout12_nLP7_dR0.2.txt')
   #Solve
   m.solve(output_dir='./test_outputs_1')
   
The user can test this code out by running the jupyter notebook `here <https://github.com/arthijayaraman-lab/crease_ga/blob/master/tutorial/workflow-simplified.ipynb>`_. 


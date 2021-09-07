Contribute by developing your own plugins
=========================================

Crease_ga allows users to develop and contribute their own :ref:`"shapes" <section-shape>` through external plugins, and crease_ga is equipped to allow `external plugins to be discoverable through package metadata <https://packaging.python.org/guides/creating-and-discovering-plugins/>`_. See `here <https://github.com/zijiewu3/cga_bbvesicle>`_ for an example plugin (cga_bbvesicle) that implements the analytical model for lipid vesicles in the `Brzustowicz and Brunger (2005) paper <https://doi.org/10.1107/S0021889804029206>`_.

How should I design my plugins?
-------------------------------

The plugins should be classes that are set up exactly as one of `the builtin shapes <https://github.com/arthijayaraman-lab/crease_ga/tree/master/crease_ga/shapes>`_. In the source code of the package, a :code:`scatterer_generator` class should exist and contain a :code:`converttoIQ` method that takes as input a list of q-values (:code:`qrange`) and an individual (:code:`param`), and returns a scattering profile I(q) evaluated at these q-values. Other helper methods can be added as you like, but :code:`converttoIQ` is the bare minimum and the only method that will be directly called by crease_ga. It is also advised that for each shape-specific descriptor and input parameter, the definition, default value (for shape-specific parameters), and default min/max bounds (for input parameter) should be explained as docstrings for the :code:`scatterer_generator` class.

What do I need to do to make my plugin discoverable?
----------------------------------------------------

You will simply need to specify the entry point for your :code:`scatterer_generator` class in your :code:`setup.py`. 

A simplified example of :code:`setup.py` is as following:

.. code:: ipython3
    
    from setuptools import setup

    setup(
        name="name-of-package",
        install_requires="crease_ga",
        entry_points={"crease_ga.plugins":["registered-name-of-plugin=name-of-package.scatterer_generator:scatterer_generator"]},
        py_modules=["name-of-package"],
        )

This will allow the plugins to be registered with crease_ga by the name :code:`registered-name-of-plugin`. If both crease_ga and your plugin are installed properly (in the order of crease_ga first, then your plugin), the :code:`scatterer_generator` class of your plugin can be loaded with the following code:

.. code:: ipython3
    
    import crease_ga.plugins as plugins
    a_name_you_prefer = plugins['registered-name-of-plugin'].load()
    new_scatterer_generator = a_name_you_prefer()

More commonly, instead of explicitly loading the scatterer_generator, all you will need is to load the shape from your plugin to `crease_ga.Model` by using `crease_ga.Model.load_shape(shape='registered-name-of-plugin')`.

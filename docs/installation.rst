Installation
============

To install this package on a linux or macOS machine, follow these steps:


#. 
   We will be using Anaconda to handle the python environment. First, we need to `install anaconda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html>`_ on our machine. You can also consider installing `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ for faster installation and lower storage consumption *only if you are an advanced user*. We recommend installing the full version of anaconda if you have limited experience with python and/or computer programming.

#. 
   At this point, you will need a bash terminal to work with. All the installation steps *after* this step will need to be accomplished from a terminal, using command line interface. 


   * 
     If you are on a linux or MacOS machine

      You can directly launch a terminal.

   * 
     If you are on a Windows machine

      If you have installed the full version of Anaconda/Miniconda in Step 1, the most straightforward way to launch a terminal will be using the *Anaconda prompt* that comes with the conda installation. You should be able to find the *Anaconda prompt* from your Start menu. You can also consider installing `Windows Subsytem for Linux (WSL) <https://ubuntu.com/wsl>`_.

#. 
   Download the package. 


   * If you have git installed, this can be done using the following command from a terminal:
     
     .. code-block::

        git clone https://github.com/arthijayaraman-lab/crease_ga

   * You can also directly download the package as a ZIP file from `our github webpage <https://github.com/arthijayaraman-lab/crease_ga>`_ by following the guidance `here <https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository-from-github/cloning-a-repository>`_. If you are following this route, you will need to unzip the package to your desired location after the download.

#. 
   Create a new conda environment and install the package along with all the dependencies. 


   * *Navigate into the root directory of the cloned package*. If you are using the anaconda prompt, you can look up `common windows command line prompts <http://www.cs.columbia.edu/~sedwards/classes/2015/1102-fall/Command%20Prompt%20Cheatsheet.pdf>`_. If you are using a unix-based shell (linux, macOS, or WSL subsystem), you can look up `common commands for Unix <http://www.mathcs.emory.edu/~valerie/courses/fall10/155/resources/unix_cheatsheet.html>`_. Either case, all you will need would probably be *displaying list of files and directories in the current folder*\  (\ ``dir`` for windows, ``ls`` for unix), and *moving to different directories*\  (\ ``cd [new_directory_path]`` for both windows and unix). You should end up at a directory called ``crease_ga``\  and be able to see files named ``setup.py``\ , ``environment.yml`` and ``README.md`` (among others) in the directory.
   * create a fresh conda environment with creaes_ga package and its dependencies installed using
     
     .. code-block::

        conda env create -f environment.yml

   * Activate the environment with
   
     .. code-block::

        conda activate crease_ga

To check if the package and its dependencies have been properly installed, you can try running

   .. code-block::

      python3

   to launch python, and then in the resulting python command line, run

   .. code-block::
      
      import crease_ga
      crease_ga.__version__
   
   If everything is properly installed, you should see the current version number of the package printed. In that case, you are all set to use the ``crease_ga`` package! Remember to activate the proper environment every time by with ``conda activate crease_ga``. 
   
   You can run the Jupyter notebook tutorials with the command

   .. code-block::

       jupyter notebook

* 
  **NOTE**\ : if you intend to run this on a supercomputing cluster, you will need to follow the steps to create a python environment on the corresponding cluster.

Try the package without installation
____________________________________

* 
  If you would like to first try our package by running our tutorial, you can directly launch a docker image of our environment, and access and interact with our jupyter notebook tutorial from a web browser *without performing any installation steps above*. To do this,  click this badge:

  .. image:: https://mybinder.org/badge_logo.svg
     :target: https://mybinder.org/v2/gh/arthijayaraman-lab/crease_ga/master
     :alt: Binder


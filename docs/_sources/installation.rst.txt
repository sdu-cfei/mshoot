.. _Installation:

============
Installation
============

Curently, the only way to install `mshoot` is by cloning its repository and installing using pip
from the local directory::

    git clone https://github.com/sdu-cfei/mshoot.git mshoot
    cd mshoot
    pip install .

Dependencies: `numpy`, `pandas`, `scipy`, `matplotlib`, `scikit-learn`, `pyfmi`.

Advised procedure
=================

Due to reliance on `pyfmi`, it is advised to use conda. This package is a bit tricky
to be installed via pip. In the below procedure it is assumed that you already have git_.

.. _git: https://git-scm.com/

1. Create a new environment: ``conda create --name mshoot python=3.6``
2. Activate the environment: ``conda activate mshoot``
3. Install dependencies from your default channel: ``conda install numpy pandas scipy matplotlib scikit-learn``
4. Install `pyfmi` from conda-forge: ``conda install -c conda-forge pyfmi``
5. Clone `mshoot` repository: ``git clone https://github.com/sdu-cfei/mshoot.git mshoot``
6. Switch directory: ``cd mshoot``
7. Install `mshoot`: ``pip install -e .``
8. Test the installation: ``python run_test.py``

Afterwards, each time you want to update mshoot to the newest version from github, just run ``git pull`` from the repository directory.
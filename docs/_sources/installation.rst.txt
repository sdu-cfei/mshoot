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
to be installed via pip.

1. Create a new environment: ``conda create --name mshoot python=3.6``
2. Activate the environment: ``conda activate mshoot``
3. Install dependencies: ``conda install numpy pandas scipy matplotlib scikit-learn pyfmi``
4. Clone `mshoot` repository: ``git clone https://github.com/sdu-cfei/mshoot.git mshoot``
5. Switch directory: ``cd mshoot``
6. Install `mshoot`: ``pip install .``
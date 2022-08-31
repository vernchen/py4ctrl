# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 09:32:18 2022

@author: Vern Chen

                  __    __    __    __
                 /  \  /  \  /  \  /  \
                /    \/    \/    \/    \
███████████████/  /██/  /██/  /██/  /████████████████████████
              /  / \   / \   / \   / \  \____
             /  /   \_/   \_/   \_/   \    o \__,
            / _/                       \_____/  `
            |/
        ███╗   ███╗ █████╗ ███╗   ███╗██████╗  █████╗
        ████╗ ████║██╔══██╗████╗ ████║██╔══██╗██╔══██╗
        ██╔████╔██║███████║██╔████╔██║██████╔╝███████║
        ██║╚██╔╝██║██╔══██║██║╚██╔╝██║██╔══██╗██╔══██║
        ██║ ╚═╝ ██║██║  ██║██║ ╚═╝ ██║██████╔╝██║  ██║
        ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═════╝ ╚═╝  ╚═╝

        mamba (0.22.1) supported by @QuantStack

        GitHub:  https://github.com/mamba-org/mamba
        Twitter: https://twitter.com/QuantStack

█████████████████████████████████████████████████████████████

1. Mambaforge can be got from:
    https://github.com/conda-forge/miniforge
    
2. After install, run below commands:
	mamba install control
	mamba install numpy
	mamba install matplotlib
	mamba install scipy
	mamba install sympy
	mamba install slycot
	mamba install spyder
	mamba install jupyterlab


Miniforge allows you to install the conda package manager with the following 
features pre-configured:

    > conda-forge set as the default (and only) channel.
        * Packages in the base environment are obtained from the conda-forge channel.
    > Optional support for PyPy in place of standard Python interpreter (aka "CPython").
    > Optional support for Mamba in place of Conda.
    > An emphasis on supporting various CPU architectures (x86_64, ppc64le, and 
                                                           aarch64 including Apple M1).

It can be compared to the Miniconda installer.

Mamba is a reimplementation of the conda package manager in C++.

    > parallel downloading of repository data and package files using multi-threading
    > libsolv for much faster dependency solving, a state of the art library used 
      in the RPM package manager of Red Hat, Fedora and OpenSUSE
    > core parts of mamba are implemented in C++ for maximum efficiency
At the same time, mamba utilize the same command line parser, package installation 
and deinstallation code and transaction verification routines as conda to stay as 
compatible as possible.

Mamba is part of a bigger ecosystem to make scientific packaging more sustainable. 
You can read our announcement blog post. The ecosystem also consists of quetz, 
an open source conda package server and boa, a fast conda package builder.

"""

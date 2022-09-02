# py4ctrl -- python for control of power electronics.


##### Target of this project is to help understanding theory of power electronics, including: Modeling, Simulation and so on.

Currently Basic buck, boost, and buck–boost converters with current programmed control are summarized, as well as Control-to-output and line-to-output transfer functions for more accurate model.

Resonant Converters design analysis, take LLC as example, simulated voltage gain function, Load-Dependent Properties of Resonant Converters and tank transfer function. 

Also included theory of automatic control system coded in python, such as Fourier Series and Fourier Transform, basic types of compensators / filters, relationship between s-domain and z-domain.

##### Recommend to run this project under mamba enviroment, with spyder installed:

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
	
##### For example, below plot output by TankNet module shows ZVS/ZCS boundary of LLC resonant tank:

![Figure 2022-09-02 192552](https://user-images.githubusercontent.com/26539320/188131845-78d4e9d3-765d-4813-9360-441ccba362e5.png)

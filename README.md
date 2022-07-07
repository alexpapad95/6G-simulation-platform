# simulation_platform
This simulation platform is presented in An Open Platform for Simulating the Physical Layer of 6G Communication Systems with Multiple Intelligent Surfaces.In this paper, we develop an open simulation platform, aimed at modeling the physical-layer electromagnetic coupling and propagation between RIS pairs. We present the platform by initially designing a basic unit cell, and then proceeding to progressively model and simulate multiple and larger RISs. The platform can be used for producing verifiable stochastic models for wireless communication in multi- RIS deployments, and its code is freely available to the public.

# Description
The simulation platform consists of two components; fine_tuning.m & Ris_pairs_simulation.m. In both parts, we use an Gaussian excitation. In boundary conditions, the MUR ABSORBING has been preferred against the PML. This selection is a fair trade-off between accuarcy and short run-time of simulation. The user is able to implement PML via the dedicated variable.

## Fine tuning
The user is able to re-adjust:  
    
    1.The central frequency and the bandwidth.

    2.The dimensions of unit cell and the thickness of S-SRRs metal.

    3.The properties (electric permittivity, εr , and tangent loss,tand) and the dimensions of substrate (width, length and thickness).

    4.The distance between the pair of S-SRRs in any axis. The distance should be kept stable in the fine tuning procedure.

The meshing is created automatically. The exported data are the following:

    1.The internal impedance for all the ports, either active or passive.

    2.The power incoming, reflected and accepted for all the ports either active or passive. 

    3.The s-parameters in the complex format for all the ports either active or passive.


## Ris_pairs_simulation

The code Ris_pairs_simulation is simulated two identical RISs.
Each RIS consists of identical S-SRRs that are examined in fine_tuning.m file.
The position of the RISs could be changed through the RISi_k variables. 
The dimensions of the RISs could be changed through the dim_meta surfaces.
Only odd numbers must be selected. The maximum value that the meshing structure is able to support is 13.
For values lower than 13,the meshing structure and the exported data are adopted automatically. We utilize the function tooclose in all the axis, reducing the simulation run-time. 
For greater dimensions, the meshing structure must be adjusted. One way is the reduction of both max_res from lamda/40 to lamda/20 
and coarseResolution from lamda/20 to lamda/10.
Another apporach is the adjustment of the tooclose function, increasing the respective criterion in "find" function. 

The exported data are the following:

    1.The internal impedance for all the ports, either active or passive.

    2.The power incoming, reflected and accepted for all the ports either active or passive. 

    3.The s-parameters in the complex format for all the ports either active or passive.

    4.The frequency resonation of the active ports. This matrix includes the position of the frequency in the respective array.

    5.The frequency in which is active ports appears the stornger coupling among any passive ports. This matrix includes the position of the frequency in the respective array.


# Steps
The user must install : MATLAB/Octave & openEMS (https://openems.de/start/)

The steps are the following:

    1.Open the MATLAB/Octave environent

    2.Open openEMS via the addpath command (it is described in the tutorial)

    3.Open the fine_tuning.m file and design the unit cell in the selected resonating frequency.

    4.Pass in the Ris_pairs_simulation the dimensions of the unit cell. In the beginning, work on the 3x3 format in order the dimensions of patches to be resulted.

    5.Upgrade the dimension of RISs keeping all the other variables constant.

The run-time is based on the selected dimension of the RISs and the computational resources.

#  Simulation Platform
This simulation platform is presented in "An Open Platform for Simulating the Physical Layer of 6G Communication Systems with Multiple Intelligent Surfaces". In this paper, we develop an open simulation platform, aimed at modeling the physical-layer electromagnetic coupling and propagation between RIS pairs. We present the platform by initially designing a basic unit cell, and then proceeding to progressively model and simulate multiple and larger RISs. The platform can be used for producing verifiable stochastic models for wireless communication in multi- RIS deployments, and its code is freely available to the public.

# Description
The simulation platform consists of two components; unit_cell_definition.m & RIS_pairs_simulation.m. Initially, unit_cell_definition.m is utilized for the design of the unit cell. The unit cell consists of the substrate, the groundplane and the S-SRR. The properties of the substrate are user-defined while the groundplane and the S-SRR are layers of metal. In the gap of S-SRR's outer ring, a lumped port is positioned. This port could radiate or not. The radiation is a Gaussian excitation
In both parts, we use an Gaussian excitation. In boundary conditions, the MUR ABSORBING has been preferred against the PML. This selection is a fair trade-off between accuarcy and short run-time of simulation. The user is able to implement PML via the dedicated variable.
In the next step, the RIS_pairs_simulation.m creates two identical metasurfaces composed by periodically positioned, identical, unit cells. The S-SRRs are connected with the adjacent ones via load patches. The distance between the RIS pairs is, also, user-defined.

## Unit cell definition
The user is able to re-adjust:  

   - The central frequency and the bandwidth of the Gaussian excitation. 
               
   - The dimensions of S-SRRs.
   
   - The properties (electric permittivity, εr, and tangent loss,tand) and the dimensions of substrate (width, length and thickness).
   
   - The distance between the pair of unit cells. By default, the distance is 10 mm.


## RIS_pairs_simulation

The first action in the RIS_pairs_simulation.m is the definition of the S-SRRs dimensions as they have resulted from the unit_cell_definition.m procedure. Subsequently, the user is able to re-adjust:

   - The dimension of the RIS pairs. The dimension displays the number of the unit cells per row and column.
   
   - The distance between the RIS units.
   
   - The distance between the S-SRRs.
   
   - The width and the thickness of the load patches.


The simulation platform can support only odd values for the dimension of the RIS (dim_meta). For values lower than 15, the meshing structure is adopted automatically. We utilize the function tooclose in all the axis, reducing the simulation run-time. For greater dimensions, the meshing structure must be adjusted. One way is the reduction of both max_res and coarseResolution from lamda/40 and lamda/20 to lamda/20 and lamda/10, correspondingly.
Another apporach is the adjustment of the "tooclose" function, increasing the respective criterion in "find".  

The exported data are the following:

The simulation is calculated the following outputs:

- The feed point impedance of the active ports.
- The incoming, reflected and accepted (subtraction of incoming and reflected) power in both active and passive ports.
- The reflection coefficients of the active ports and the transmission coefficients of them with the passive ones (s-parameters).
- The resonating frequency in which the reflection coefficients of the active ports are minimized.
- The values of frequency in which the active ports appear the maximized transmission coefficients with the passive ones.



# Steps
The user must install : MATLAB/Octave & openEMS (https://openems.de/start/)

The steps are the following:

   - Open the MATLAB/Octave environent

   - Open openEMS via the addpath command (it is described in the tutorial)

   - Open the unit_cell_definition.m file and design the unit cell in the selected resonating frequency.

   - Pass in the RIS_pairs_simulation the dimensions of the unit cell. In the beginning, work on the 3x3 format in order the dimensions of load patches to be resulted.

   - Upgrade the dimension of RISs keeping all the other variables constant.

The run-time is based on the selected dimension of the RISs and the computational resources.

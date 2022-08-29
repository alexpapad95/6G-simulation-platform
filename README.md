#  6G Simulation Platform: From the Physical-layer to Networking with RIS

In this paper, we develop an open simulation platform, aimed at modeling the physical-layer electromagnetic coupling and propagation between RIS pairs. We present the platform by initially designing a customizable physical layer, comprising a basic unit cell, and then proceeding to progressively model and simulate multiple and larger RISs. Then, we describe how the platform can be used for producing verifiable stochastic models for wireless communication in multi- RIS deployments, facilitating realistic network-layer studies of RIS-enabled 6G networks. 

The source code is freely available.

The platform is presented in the following publication(s):

   - Alexandros Papadopoulos, Antonios Lalas, Konstantinos Votis, Dimitrios Tyrovolas, George Karagiannidis, Sotiris Ioannidis, Christos Liaskos."An Open Platform for Simulating the Physical Layer of 6G Communication Systems with Multiple Intelligent Surfaces", In Proceedings of CNSM 2022,18th International Conference on Network and Service Management.

# Installation Instructions
The user must install : MATLAB/Octave & openEMS (https://openems.de/start/)

The steps are the following:

   - Open the MATLAB/Octave environment

   - The openEMS tool is activated through the followind command:
   ``` 
   addpath ('C:\openEMS\matlab');
   ```
   - Open the unit_cell_definition.m file and design the unit cell in the selected resonating frequency.

   - Pass in the RIS_pairs_simulation.m the dimensions of the unit cell. In the beginning, work on the 3x3 format in order the dimensions of load patches to be resulted.

   - Upgrade the dimension of RISs keeping all the other variables constant.

The run-time is based on the selected dimension of the RISs and the computational resources. In our experiments, we utilized a PC with 64-bit Windows 10, installed RAM 32GB and processor Intel(R) Core(TM) i7-10750H CPU @ 2.60 GHz. For the given specifications of the device, the run-times are presented in the following matrix:
 
 <div align="center">

| Dimensions of RIS   |      Run-time of code       
|----------|:-------------:|
| 1x1 unit cells| 1min |
| 3x3 metasurfaces | 2.2min |
| 11x11 metasurfaces | 1.5 hours |

</div>

The platform could be supported by multicore processing but not by GPU one. This limitation is sourced from the openEMS tool. 

# Physical Layer Simulation: General Description 
The proposed platform models the accurate electromagnetic propagation between two RIS units–whose composition, dimensions, and state are user-defined–while varying
the distance between them. Each RIS consists of a planar arrangement of elements, known as unit cells. Each element hosts a port and a lumped element offering tunable impedance within a user-specified range. We consider that all the elements of the RIS1 are active, meaning that power is generated at each of their ports, enters the system and is altered by the state of the corresponding lumped element. Additionally, all the elements of the RIS2 are passive, meaning that their ports
receive power emanated from RIS1, alter it based on the state of the local lumped element and re-emit it. Depending on its configuration, this setup emulates either: i) the RIS-RIS part of a transmitter-RIS-RIS-receiver communication, which is useful for deriving verifiable stochastic models for the inter-RIS part of a channel. ii) the transmitter-RIS or the RIS-receiver part of the communication.

The simulation platform consists of two components; unit_cell_definition.m & RIS_pairs_simulation.m. 

 <div align="center">

   ![1](https://user-images.githubusercontent.com/72256279/186215268-bfe07550-6c0b-4da2-b7f9-911d88d0451e.jpg)

 
 </div>  
 
The unit cell consists of the substrate, the groundplane and the Square-Split Resonator Ring (S-SRR). The properties of the substrate are user-defined while the groundplane and the S-SRR are layers of metal. In the gap of S-SRR's outer ring, a lumped port is positioned. This port could radiate or not. The radiation is a Gaussian excitation. 


 <div align="center">
   
   ![3](https://user-images.githubusercontent.com/72256279/186148890-ae054ed0-d126-4454-9935-8d575735a5aa.png)

</div>  

Initially, the unit cell definition is utilized. The user determines the resonating frequency, meaning the frequency of the Gaussian excitation, the properties of the substrate and the positions, consequently and the distance, of the unit cells. These parameters must be remained stable in the following procedure. After that the fine tuning stage begins. The user re-defines the dimensions of both substrate and S-SRR until the resonation to be achieved. The criteria that should be held in order the fine tuning to be done are described in the Section III, subsection C of the paper. The whole procedure is presented in Section IV/Step 1.

After that, the user passes the dimensions of the substrate and the S-SRRs in the RIS_pair_simulation. Apart from that, the dimension of the RISs and the position of them have also to be defined and remained stable. The fine tuning procedure in this phase includes the determination of the distance between the adjacent S-SRRs and the dimensions of the load patches. This process is presented in Section IV/Steps 2,3. Once it is completed, the user is able to change the position of the RISs observing their behavior based on the output data (see below).


## Unit cell definition
The user is able to re-adjust the following variables:

   ``` 
        f0=8e9; % central frequency of the Gaussian excitation
        fc=2e9  % bandwidth of the Gaussian excitation
       
       %% dimensions of S-SRRs %%
        L1=10.9; % length of the outer ring
        L2=10.36; % width of the outer ring
        G1=1.05; % width of the outer gap
        width_outer=1.05; % width outer patch
        width_iner=width_outer;
        G3=1.05; %distance between iner and outer
        G2=G1; % width of inner gap
        srr_thickness=0.15;
       
       %% properties of substrate %%
        substrate.epsR = 2.2; % electric permittivity
        tand = 0.024; % tangent loss
        substrate_srr.kappa = tand*2*pi*f0*EPS0*substrate.epsR;
        feed.R = 50; % feed resistance
       
       %% dimensions of substrate %%
        dx=1.85; 
        dz=1.85+L2-L1;
        substrate_srr.width = L2+dx;
        substrate_srr.length =L1+dz;
        substrate_srr.thickness = 4.8;

       %% S-SRRs positions and distance between them %%
        ssrr_x=0;
        ssrr_y=5;
        ssrr_z=0;
       
        ssrr2_x=0;
        ssrr2_y=-5;
        ssrr2_z=0;
   ```


## RIS_pairs_simulation

The first action in the RIS_pairs_simulation.m is the definition of the S-SRRs dimensions as they have resulted from the unit_cell_definition.m procedure. Subsequently, the user is able to re-adjust:
   ```
       %% dimensions of RISs %%
        dim_meta=3; % The dimension displays the number of the unit cells per row and column.
       
       %% distance between the RIS units in any axis %%
        RIS1_x=0;
        RIS1_y=5;
        RIS1_z=0;

        RIS2_x=0; 
        RIS2_y=-5;
        RIS2_z=0;
       
       %% distance between the S-SRRs unit cells %%
        dx=10; % x-axis distance between elements
        dz=10; % z-axis distance between elements

       %% width and thickness of the load patches %%
        width_patch=1.2;
        patch_thickness=srr_thickness;
       
   ```

## Output Data
The simulation produces the following outputs:

- The feed point impedance (Zin) of the active ports.
- The incoming (P_incoming), reflected (P_reflected) and accepted (P_accepted) (subtraction of incoming and reflected) power in both active and passive ports. The accepted power is the substration of the incoming and reflected power. Therefore, it may be negative for the passive ports.
- The reflection coefficients of the active ports and the transmission coefficients of them with the passive ones (s-parameters).
- The resonating frequency (f_res) in which the reflection coefficients of the active ports are minimized.
- The values of frequency in which the active ports appear the maximized transmission coefficients with the passive ones (f_coupling).

There are two entities that determine the size of the output data; the dimension of RISs (dim_meta) and the frequency interval in which the simulation platform works. We define the interval as (f0-fc,f0+fc) and we separate it in 1001 parts. The dedicated command is presented below. By all means, the user is able to increase or decrease the number of parts for higher accuracy or for lower memory demand, correspondingly.

   ```
      freq = linspace(f0-fc,f0+fc,1001);
   ```
In the following matrix, we present the dimension of each output using two parameters; dim_meta and interval_parts.

 <div align="center">

| Output data   | Dimension     
|----------|:-------------:|
| Zin | 1x 2dim_meta $^2$ |
| P_incoming | 1 x 2dim_meta $^2$ |
| P_reflected | 1 x 2dim_meta $^2$ |
| P_accepted | 1 x 2dim_meta $^2$ |
| s-parameters | 2dim_meta $^2$ x 2dim_meta $^2$ x interval_parts |
| f_res | 1 x dim_meta $^2$ |
| f_coupling | 2dim_meta $^2$ x 2dim_meta $^2$ |

</div>

The output data can be saved as MATLAB files or as CSV files using the following commands:

   ```
     writecell(f_res,'f_res.csv')
     writecell(f_coupling,'f_coupling.csv')
     writecell(s,'s_parameters.csv')
     writecell(Pincoming,'P_incoming.csv')
     writecell(Preflected,'P_reflected.csv')
     writecell(Paccepted,'P_accepted.csv')
     writecell(P_in,'P_in.csv')
   ```
Some indicative output data for dim_meta=3 and interval_parts=1000 are the following:
[P_in.csv](https://github.com/alexpapad95/6G_simulation_platform/files/9400473/P_in.csv), 
[P_incoming.csv](https://github.com/alexpapad95/6G_simulation_platform/files/9400474/P_incoming.csv), 
[P_reflected.csv](https://github.com/alexpapad95/6G_simulation_platform/files/9400476/P_reflected.csv), 
[P_accepted.csv](https://github.com/alexpapad95/6G_simulation_platform/files/9400478/P_accepted.csv), 
[s_parameters.csv](https://github.com/alexpapad95/6G_simulation_platform/files/9400466/s_parameters.csv), 
[f_coupling.csv](https://github.com/alexpapad95/6G_simulation_platform/files/9400468/f_coupling.csv), 
[f_res.csv](https://github.com/alexpapad95/6G_simulation_platform/files/9400470/f_res.csv). 

## Meshing in openEMS
The simulation platform can support only odd values for the dimension of the RIS (dim_meta). For values lower than 15, the meshing structure is created automatically. We utilize the function "tooclose" in all the axis, reducing the simulation run-time. For larger dimensions, the meshing structure must be adjusted via 2 strategies that could be utilized simultaneously. 

### Strategy 1
The first action is the decrease of the resolution values. The initial values are:
 ```
        max_res = c0 / (f0 + fc) / sqrt(substrate_srr.epsR) / unit /40; 
        coarseResolution = c0/(f0 + fc) / unit / 20;
   ```
The value for larger dimensions of RISs could be:
 ```
        max_res = c0 / (f0 + fc) / sqrt(substrate_srr.epsR) / unit /20; 
        coarseResolution = c0/(f0 + fc) / unit / 10;
   ```
### Strategy 2
The "tooclose" function intervenes to avoid a long-time simulation in the case that the distance between adjacent meshing lines becomes less than a pre-defined
threshold in each axis. By default, this threshold is equal to the half of the maximum resolution. 
 ```
tooclose = find (diff(mesh.y) < max_res/2);
if ~isempty(tooclose)
    mesh.y(tooclose) = (mesh.y(tooclose) + mesh.y(tooclose+1))/2;
    mesh.y(tooclose + 1) = [];
end

tooclose = find (diff(mesh.x) < max_res/2);
if ~isempty(tooclose)
    mesh.x(tooclose) = (mesh.x(tooclose) + mesh.x(tooclose+1))/2;
    mesh.x(tooclose + 1) = [];
end

tooclose = find (diff(mesh.z) < max_res/2);
if ~isempty(tooclose)
    mesh.z(tooclose) = (mesh.z(tooclose) + mesh.z(tooclose+1))/2;
    mesh.z(tooclose + 1) = [];
end
   ```
The second strategy is the increase of the function's threshold in order it to intervene more often. The code could be configured as:
 ```
tooclose = find (diff(mesh.y) < max_res);
if ~isempty(tooclose)
    mesh.y(tooclose) = (mesh.y(tooclose) + mesh.y(tooclose+1))/2;
    mesh.y(tooclose + 1) = [];
end

tooclose = find (diff(mesh.x) < max_res);
if ~isempty(tooclose)
    mesh.x(tooclose) = (mesh.x(tooclose) + mesh.x(tooclose+1))/2;
    mesh.x(tooclose + 1) = [];
end

tooclose = find (diff(mesh.z) < max_res);
if ~isempty(tooclose)
    mesh.z(tooclose) = (mesh.z(tooclose) + mesh.z(tooclose+1))/2;
    mesh.z(tooclose + 1) = [];
end
   ```


# Usage in wireless channel modeling and networking 
The authors intend to update the repository for a span of three studies on the subject of RIS-enabled 6G communications as follows:
1. Physical layer, covering the RIS pair design and EM simulation, up to the derivation of S-parameters. (Present paper)
2. Macroscopic channel modeling (path loss and fading) for inter-RIS communications, based on multiple runs of the present platform, on various settings (indoors, outdoors, mobility, etc.) and using stochastic modeling.
3. Network layer, employing the statistical models of #2 for the real-time simulation of a communication network, studying the effect of RIS on MAC and network-layer protocols.
As the project progresses, additional sections and source code will be added for steps #2 and #3.
Here we present an outline of this process. 
First, consider the S-parameters exported by the platform for a given RIS pair configuration (i.e., specific EM design and distance). Each specific parameter describes the amount of power (and phase) entering each element of the receiving RIS, from each element of the transmitting RIS. Notice that the employed approach, i.e., using the S-parameters as the basic data export unit of the platform, provides the benefit of easily cascading multiple, successive RIS units, and deducing the received signal at the end-point(s) of interest. 
Second, the macroscopic channel modeling seeks to extract a stochastic channel model for each setting of interest, in order to avoid the need for frequent and time-consuming platform runs.The channel model considers various random RIS placements (varying distance, orientation, position in a real setup).For each randomized case, multiple platform runs are executed and stochastic models are extracted to describe the path loss and multipath fading phenomena per case.
(Notice that this step can be skipped, and instead invoke runs in the platform as required for a specific setup.However, the existence of a validated stochastic model adds versatility to the simulation process).
Finally, the network layer abstracts the physics of RIS and treats the system as a graph. The graph nodes are the RIS units and the user devices. A link is inserted in the graph for every node pair in connectivity (e.g., line of sight). Each link is assigned a macroscopic wireless channel model from step #2, based on the conditions around the node pair (e.g., mobility, distance, obstacles around them, etc.). In tandem with an encoding and modulation of interest, the data error model and throughput can be deduced per link. From this point and on, network-layer analysis can be performed, which can model other network infrastructure (e.g., routers) as additional nodes in the graph. Notice that the network layer can model the control latency between the RIS units and the central (e.g., SDN) controller that continuously reconfigures them in order to constantly adapt to the user objectives and the overall system state.



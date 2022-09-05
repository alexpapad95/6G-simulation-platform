#  6G Simulation Platform: Examples
In this folder, some examples of the platform utilization are included. First, we create unit cells and RIS pair at 5 and 8GHz with the usage of the Rogers RTDuroid 5880 substrate. The respective codes are:
-  unit_cell_definition_5GHz_Rogers_RTDuroid 5880.m
-  unit_cell_definition_8GHz_Rogers_RTDuroid 5880.m
-  RIS_pairs_simulation_5GHz_Rogers_RTDuroid 5880.m
-  RIS_pairs_simulation_8GHz_Rogers_RTDuroid 5880.m

Additionally, we present some examples in which lumped elements are included in order the local tunability to be feasible. The dedicated codes are:
- RIS_pairs_simulation_8GHz_tunable_ports.m
- RIS_pairs_simulation_8GHz_tunable_load_patches.m


## Resonating in different frequencies with different substrate types 

We have created RIS pairs at 5GHz with the usage of Rogers RTDuroid 5880 substrate. Firstly, we define the resonating frequency at 5GHz and the properties of substrate. We complete the fine tuning procedure in unit_cell_definition_5GHz_Rogers RTDuroid 5880.m as it is described in the 6G-simulation-platform. After, that we pass the dimensions of the S-SRR in the RIS_pairs_simulation_5GHz_Rogers RTDuroid 5880.m and we work with dim_meta=3 in order the final setup to be determined. 
The same procedure is repeated for 8GHz in the codes unit_cell_definition_8GHz_Rogers RTDuroid 5880.m and RIS_pairs_simulation_8GHz_Rogers RTDuroid 5880.m

The user is able to select any substrate types.The most widespread types of them that are described in the research literature are the following:

 <div align="center">

| Name of material   |      dielectric constant   |  tangent loss
|----------|:-------------:|:-------------:
| Rogers RO4350B | 3.66  |   0.0037 |
| Rogers RTDuroid 5880 | 2.2 | 0.0009 |
| FR-4 | 4.3 | 0.025 |
| Lossy Polymer | 2.9  | 0.00069 |

</div>

## Local tunability of the RISs


In the beginning steps of metasurfaces research, once the unit cell is constructed, both its function and the operation frequency are completely fixed. For instance, an absorber operates at a specific, pre-defined, frequency where the input impedance is matched to the free space. 

However, because of the structural makeup of the unit cells, re-design procedures are required if we want to reconfigure the operating frequency and/or its functionality. The main functionalities are beam-focusing, beam-splitting, change of polarization and filtering. 

By including "tuning" capacity in the unit cells, the properties of the metamaterials and metasurfaces can be changed. Then, by altering the stimulus, their electromagnetic wave behavior can be modified externally. The tuning can also be managed by a computer program. According to this perspective, these programmable metasurfaces offer greater possibilities for achieving dynamical wave applications without the need of re-fabrication.

The tunability of metasurfaces' unit cells could be either global or local. In the global tunability, all the states of the unit cells are reconfigured with the same way. This method is far easier than the local one but lacks of the functionalities that could be achieved. 

Simulating a local tunability mechanism, we have introduced in the previous codes lumped elements that could be configured with different way in two points; in the top of the lumped ports that radiate and the at middle of the load patches that connect the adjacent S-SRRs. The common utilization of both methods could easily be done. 

We use values for reactance from 0 to 5 Ohms and for capacitance from 1 to 5 pF, corresponding to the parameters of commercially available chip varactors for the target frequency bands.

We work on the structure that has been resulted from the usage of the Rogers RTDuroid 5880 substrate at 8GHz for dim_meta=3. The depicted structure is created with the usage of RIS_pairs_simulation_8GHz_tunable_ports.m:

 <div align="center">
   
   ![ports](https://user-images.githubusercontent.com/72256279/188427806-34c14a1f-e0ac-48d2-9406-f83145f870b3.PNG)


</div>  

The user is able to reconfigure the values of reactance and capacitance in the ports through the dedicated commands:

``` 
lumped_c = randi([10 50], 2,dim_meta^2)/10;
lumped_c=lumped_c*1e-12;
lumped_r= randi([0 50], 2,dim_meta^2)/10;
   ```

In the same format, the usage of RIS_pairs_simulation_8GHz_tunable_load_patches.m introduce lumped elements at the middle of the patches. The structure for dim_meta=3 is the following:

 <div align="center">
   
   ![load_patches_1](https://user-images.githubusercontent.com/72256279/188427829-596cb144-a1aa-4e15-9db9-1b01b2bdc251.PNG)
   ![load_patches_2](https://user-images.githubusercontent.com/72256279/188427848-712c95f4-b468-44a7-b350-954804cf1219.PNG)
</div>  

The user is able to reconfigure the values of reactance and capacitance in the load patches through the dedicated commands:

``` 
lumped_c = randi([10 50], 2,4*dim_meta*(dim_meta-1))/10;
lumped_c=lumped_c*1e-12;
lumped_r= randi([0 50], 2,4*dim_meta*(dim_meta-1))/10;
   ```



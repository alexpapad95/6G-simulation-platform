#  6G Simulation Platform: Examples
In this folder, some examples of the platform utilization are included. 

## Resonating in different frequencies with plenty sustrate properties 

We have created RIS pairs at 5GHz and 8GHz at Rogers RTDuroid 5880 substrate. The respective codes are unit_cell_definition_5GHz_Rogers RTDuroid 5880.m and unit_cell_definition_8GHz_Rogers RTDuroid 5880.m, correspondingly, for the design of the unit cell and the required fine tuning procedure. The RIS_pairs_simulation_5GHz_Rogers RTDuroid 5880.m and RIS_pairs_simulation_8GHz_Rogers RTDuroid 5880.m create the respective pairs of the RIS.
The most widespread types of substrate that could be utilized are the following:

 <div align="center">

| Name of material   |      dielectric constant   |  tangent loss
|----------|:-------------:|:-------------:
| Rogers RO4350B | 3.66  |   0.0037 |
| Rogers RTDuroid 5880 | 2.2 | 0.0009 |
| FR-4 | 4.3 | 0.025 |
| Lossy Polymer | 2.9  | 0.00069 |

</div>

## Local tunability of the RISs

We have introduced in the previous codes lumped elements in two points; in the top of the lumped ports that radiate and the at middle of the load patches that connect the adjacent S-SRRs. The common utilization of both methods could easily be done. 
We use values for reactance from 0 to 5 Ohms and for capacitance from 1 to 5 pF, corresponding to the parameters of commercially available chip varactors for the target frequency bands.

The structure that is created with the usage of RIS_pairs_simulation_5GHz_tunable_ports.m and RIS_pairs_simulation_8GHz_tunable_ports.m, for dim_meta=3, is the following:

 <div align="center">
   
   ![ports](https://user-images.githubusercontent.com/72256279/188427806-34c14a1f-e0ac-48d2-9406-f83145f870b3.PNG)


</div>  

The user is able to reconfigure the values of reactance and capacitance in the ports through the dedicated commands:


``` 
    lumped_c = randi([10 50], 2,dim_meta^2)/10;
    lumped_c=lumped_c*1e-12;
    lumped_r= randi([0 50], 2,dim_meta^2)/10;
   ```


The structure that is created with the usage of RIS_pairs_simulation_5GHz_tunable_load_patches.m and RIS_pairs_simulation_8GHz_tunable_load_patches.m, for dim_meta=3, is the following:

 <div align="center">
   
   ![load_patches_1](https://user-images.githubusercontent.com/72256279/188427829-596cb144-a1aa-4e15-9db9-1b01b2bdc251.PNG)


</div>  

 <div align="center">
   
![load_patches_2](https://user-images.githubusercontent.com/72256279/188427848-712c95f4-b468-44a7-b350-954804cf1219.PNG)

</div>  

The user is able to reconfigure the values of reactance and capacitance in the load patches through the dedicated commands:


``` 
    lumped_c = randi([10 50], 2,4*dim_meta*(dim_meta-1))/10;
    lumped_c=lumped_c*1e-12;
    lumped_r= randi([0 50], 2,4*dim_meta*(dim_meta-1))/10;
   ```



# Customized RMR Solver 

This repository includes a customized version of the **Rapid Muscle Redundancy (RMR) solver** integrated with different formulations of **glenohumeral stability**.  

## Overview  

The original RMR solver, modeled glenohumeral stability as an **inequality constraint** that keeps the **glenohumeral contact force direction** within the stability border of the glenoid cavity, approximated as a **circle**.  

The original RMR solver is available at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8360269.svg)](https://doi.org/10.5281/zenodo.8360269)

In this customized version, we proposed different formulations to model glenohumeral stability using:  

- **Inequality constraints**, where the direction of the glenohumeral contact force was constrained within different stability border approximations:  
  - A **point**  
  - An **ellipse**  
  - A **polynomial fit** of empirical values obtained from **concavity compression tests**  
- **Penalty terms** in the objective function as:
  - **Conditional penalty**: that triggers if the direction of the contact force exceeds a specified circular stability border
  - **Planar penalty**: increases as the direction of the glenohumeral contact force points further away from the glenoid cavity border. Here, the glenoid cavity is considered as a flat plane. This plane is the plane that contains the glenoid cavity border.
  -  **Curve penalty**: increases as the direction of the glenohumeral contact force points further away from the glenoid cavity border. Here, the curvature of the glenoid cavity is approximated as spherical.

 ## Musculoskeletal Model
 All the analysis was performed on the thoracoscapular shoulder model **DOI:** [10.3389/fnbot.2019.0009](https://doi.org/10.3389/fnbot.2019.00090)

## Validation  

A dataset of upper body kinematics and glenohumeral contact forces, obtained from an instrumented shoulder prosthesis, was used to validate the estimated magnitude and direction of the glenohumeral contact force across the different stability formulations.
Dataset available at **DOI:** [10.4121/86db1d7d-13d9-4631-9c6b-1e3134a1ab38](https://doi.org/10.4121/86db1d7d-13d9-4631-9c6b-1e3134a1ab38)

## Requirements
The code in this repository was implemented and tested using:

  -  **MATLAB 2023a**: Optimization and Signal Processing Toolboxes are required. **Note:** using earlier version of MATLAB may give some error. For example, the code was tested on MATLAB R2021a and some functions like  `writelines()` in `main_analyse_dataset.m`, and `quiver()` in `RMR_analysis.m` were not available in this version. These two functions, however, are just for data saving and plotting, and can be commented out. Similar issues may arise with earlier MATLAB versions. 
     
  - **OpenSim 4.4**: [Download here](https://simtk.org/frs/?group_id=91), and Follow the instructions [here](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab) to set up the OpenSim MATLAB API.

  The code is not guaranteed to run on different versions of these softwares

## Running the Code
To solve muscle redundancy:
- Navigate to Customized RMR Solver and open main_analyse_dataset.m
- Define input flags to determine which stability formulation to be used
- Select the input model
- locate the trc files
- Set the saving path
- Run the file

To generate the graphs used in the paper or other customized graphs:
- Navigate to Personal_Results
- Customize/Run Generate_Graphs.m

## Issues
If you experience any problems running the code, feel free to open an issue and we will try help you!

## Publications
Please cite this work as:

- Hasan, I. M. I., Belli, I., Seth, A. & Gutierrez-Farewik, E. M. Modeling glenohumeral stability in musculoskeletal simulations: A validation study with in vivo contact forces. bioRxiv DOI: 10.1101/2025.05.08.652806 (2025).

<p align="center">
  <img src="https://github.com/user-attachments/assets/4accc270-59f7-44af-adb0-c02fddaaf384" width="800" >
</p>




## Funding
This work was funded by the Promobilia Foundation, and the Swedish Research Council, Stockholm, Sweden. 

<p align="center">
  <img src="https://github.com/user-attachments/assets/1ea03a52-8ffa-4abf-b911-5e29c71f41d0" width="800" >
</p>

## Authors
Customization of the original RMR solver and integration of the different stability formulations, included in this repository, were implemented by Ibrahim Mohammed I. Hasan during his Ph.D. at KTH MoveAbility, Department of Engineering Mechanics, KTH Royal Institute of Technology, Stockholm, Sweden. 


<p align="center">
  <img src="https://github.com/user-attachments/assets/c280173a-8689-47fa-8f33-fcb32693e0f1" width="500">
</p>

## License
- The code included in this repository is licensed under the Apache 2.0 license (see the `LICENSE_code` file).
- The thoracoscapular shoulder model is licensed under CC BY 4.0 Use Agreement terms. Find more information at **DOI:** [10.3389/fnbot.2019.0009](https://doi.org/10.3389/fnbot.2019.00090)
- The validation dataset is licensed under CC0 1.0 Universal terms. Find more information at **DOI:** [10.4121/86db1d7d-13d9-4631-9c6b-1e3134a1ab38](https://doi.org/10.4121/86db1d7d-13d9-4631-9c6b-1e3134a1ab38)



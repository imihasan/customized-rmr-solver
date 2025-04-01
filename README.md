# Customized RMR Solver 

This repository includes a customized version of the **Rapid Muscle Redundancy (RMR) solver** integrated with different formulations of **glenohumeral stability**.  

## Overview  

The original RMR solver, modeled glenohumeral stability as an **inequality constraint** that keeps the **glenohumeral contact force direction** within the stability border of the glenoid cavity, approximated as a **circle**.  

The original RMR solver is available at: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8360269.svg)](https://doi.org/10.5281/zenodo.8360269)

In this customized version, we integrated different formulations to model glenohumeral stability using:  

- **Inequality constraints**, where the direction of the glenohumeral contact force was constrained within different stability border approximations:  
  - A **point**  
  - An **ellipse**  
  - A **polynomial fit** of empirical values obtained from **concavity compression tests**  
- **Penalty terms** in the objective function,
  - **Conditional penalty** penalty: trigers if the direction of the contact force exceedes a speficied circular stability border
  - **Planar penalty**: increases as the direction of the glenhumeral contact force points further away from the glenoid cavity border. Here, the glenoid cavity is considered as a flat plane.This plane is the plane that contains the glenoid caity border.
  -  **Curve penalty**: increases as the direction of the glenhumeral contact force points further away from the glenoid cavity border. Here, the curvature of the glenoid cavity is approximated as spherical.

 ## Musculoskeletal Model
 All the analysis was performed on the thoracoscapular shoulder model **DOI:** [10.3389/fnbot.2019.0009](https://doi.org/10.3389/fnbot.2019.00090)

## Validation  

A dataset of upper body kinematics and glenohumeral contact forces, obtained from an instrumented shoulder prosthesis, was used to validate the estimated magnitude and direction of the glenohumeral contact force across the different stability formulations.
Dataset available at **DOI:** [10.4121/86db1d7d-13d9-4631-9c6b-1e3134a1ab38](https://doi.org/10.4121/86db1d7d-13d9-4631-9c6b-1e3134a1ab38)

## Requirements
The code in this repository was implemented and tested using:

  -  **MATLAB 2023a**: Optimization and Signal Processing Toolboxes are required
     
  - **OpenSim 4.4**: [Download here](https://simtk.org/frs/?group_id=91), and Follow the instructions [here](https://simtk-confluence.stanford.edu:8443/display/OpenSim/Scripting+with+Matlab) to set up the OpenSim MATLAB API.

  The code is not guaranteed to run on different versions of these softwares

## Running the Code
To solve muscle redundancy:
- Navigate to Customized RMR Solver and open main_analyse_dataset.m
- Define iput flags to determine which stability formulation to be used
- Select the input model
- locate the trc files
- Set the saving path
- Run the file

To generate the graphs used in the paper or other customized graphs:
- Navigate to Personal_Results
- Customize/Run Generate_Graphs.m

## Issues
If you experience ay problems running the code, feel free to open an issue and we will try help you!

## Publications
Please cite this work as:

- Ibrahim Mohammed I. Hasan, Italo Belli, Ajay Seth, and Elena M. Gutierrez-Farewik. "Modeling glenohumeral stability in musculoskeletal simulations: A validation study with \emph{in vivo} contact forces". Manuscript.


## Funding
This work was funded by the Promobilia Foundation, and the Swedish Research Council, Stockholm, Sweden. 

<p align="center">
  <img src="https://github.com/user-attachments/assets/1279341c-86dd-41ad-ba88-712d21e1af96" width="300" hspace="20">
  <img src="https://github.com/user-attachments/assets/6e0aeae9-4fcc-4795-b695-533111dd219e" width="300">
</p>

## Authors
Customization of the original RMR solver and integration of the different stability formulations, included in this repository, were implemented by Ibrahim Mohammed I. Hasan during his Ph.D. at KTH MoveAbility, Department of Engineering Mechanics, KTH Royal Institute of Technology, Stockholm, Sweden. 


<p align="center">
  <img src="https://github.com/user-attachments/assets/cf9c908b-f8f7-413e-9b7f-19e1e80450b1" width="500">
</p>




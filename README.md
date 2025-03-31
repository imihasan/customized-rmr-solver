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

Dataset of upper body **kinematics and glenohumeral contact force measurements** from an **instrumented shoulder prosthesis** were used to validate the estimated **magnitude** and **direction** of the glenohumeral contact force.  
Dataset available at **DOI:** [10.4121/86db1d7d-13d9-4631-9c6b-1e3134a1ab38](https://doi.org/10.4121/86db1d7d-13d9-4631-9c6b-1e3134a1ab38)

## Publications
Please cite this work as:


## Funding
This work was funded by the Promobilia Foundation, and the Swedish Research Council, Stockholm, Sweden. 

<p align="center">
  <img src="https://github.com/user-attachments/assets/1279341c-86dd-41ad-ba88-712d21e1af96" width="300" hspace="20">
  <img src="https://github.com/user-attachments/assets/6e0aeae9-4fcc-4795-b695-533111dd219e" width="300">
</p>

## Author
Customizations included in this repository were implemented by Ibrahim Mohammed I. Hasan during his Ph.D. at KTH MoveAbility, Department of Engineering Mechanics, KTH Royal Institute of Technology, Stockholm, Sweden. 


<p align="center">
  <img src="https://github.com/user-attachments/assets/cf9c908b-f8f7-413e-9b7f-19e1e80450b1" width="500">
</p>


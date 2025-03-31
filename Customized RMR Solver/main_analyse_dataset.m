% Script to run the Rapid Muscle Redundancy (RMR) solver on user-selected experiments.
% The user is prompted with the selection of the tasks to analyze.
% Within the script, it is possible to adjust the downsampling to be
% applied, and whether the analysis should include the glenohumeral
% constraint or not.
%
% Author: Italo Belli (i.belli@tudelft.nl) 2023
%
% Customized by: Ibrahim Mohammed I. Hasan (imihasan@kth.se) 2025
%
% Customization included integrating different formulations to model the
% glenohumeral stability problem. The original formulation in the rmr
% solver enforced an inequality constrain to maintain the glenohumeral
% contact froce vector within the glenoid cavity border approximated as a
% circle. The customized version included differnt approximation of the
% glenoid cavity border as a point, an ellipse, and a polynomial. Also a
% different strategy was introduced to model glebohumeral stability as a
% conditional and a continuous penlaty in the objective function. 

close all; clear; clc; beep off;

% Import the OpenSim libraries.
import org.opensim.modeling.*;

% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

% getting path to other folders in this repo
addpath(pathstr)
cd ..\
path_to_repo = pwd;
addpath(path_to_repo)
addpath(fullfile(path_to_repo, 'Code\Data Processing\'))


% Flags (Select whether to enforce constraints)
% enforcing continuity of the activations from one timestep to the next, 
% to respect first-order dynamics
dynamic_bounds = true;
% enforcing directional constraint on the glenohumeral joint force
enforce_GH_constraint = false;
% string indicating the name of the stability border to be enforced. 
% It can be: circle, ellipse, polynomial or point.
tag_con="ellipse";
% string indicating which penalty to use in the cost function: 
% "conditional", "planar", "curve", or "no_penalty"
tag_cost= "curve";
% execlude locked coordinates from matching
execlude_locked=0;
% execlude clavicle coordinates from matching
execl_clav=0;
% execlude scaula winging coordinate from matching
execl_wing=0;    

% where you have the experimental files (.trc)
trc_path = fullfile(path_to_repo, 'ExperimentalData\S8R\trc');

% where to save the results
Edit_tag='Test'; %name of the folder where results will be saved
mkdir(fullfile(path_to_repo, ['Personal_Results\S8R\' Edit_tag]));
saving_path = fullfile(path_to_repo, ['Personal_Results\S8R\' Edit_tag]);
%Write a text file to describe the changes or the conditions you have made
%in this trial
Edit_description="Description_Here";
writelines(Edit_description,fullfile(saving_path,"Description.txt"));

% Select model
modelFile_2kg = append(path_to_repo, '\OpenSim Models\OrthoModel_2kgWeight.osim');
model_2kg = Model(modelFile_2kg);

%Read a baseline model and re-assign the 2.4kg to the hand in the current
%model
model_base=Model(path_to_repo+"\OpenSim Models\BaseModel.osim");
model_2kg.updBodySet().get("hand").setMass(model_base.updBodySet().get("hand").getMass()+2.4);
model_2kg.updBodySet().get("hand").setMassCenter(Vec3(0, -0.03, 0));
model_2kg.finalizeConnections();

% Select the experimental data to be considered
dataset_considered = 'Orthoload';

[files,path] = uigetfile('*.trc', 'Select the .trc files to analyse', trc_path, 'MultiSelect','on');

if iscell(files)
    num_files = size(files, 2);
else
    num_files = 1;
end

% Set the weight for the various scapula coordinates in IK
% This is to achieve a good agreement between scapula upward rotation and
% shoulder elevation (as reported in the paper)
weight_abd = 0.0001;
weight_elev = 0.0001;
weight_up_rot = 0.0002;
weigth_wing = 0.0001;
weight_coord = [weight_abd, weight_elev, weight_up_rot, weigth_wing];

% Downsampling
time_interval = 1;
t_end=4; %set the end time

%% Run Rapid Muscle Redundancy (RMR) solver
% preallocating arrays to hold information about the solutions
optimizationStatus = [];
unfeasibility_flag = [];
tOptim = zeros(num_files,1);
result_file_RMR = {};

for trc_file_index=1:num_files
    fprintf('Running RMR on experiment %i \n', trc_file_index)
    if num_files>1
        experiment = append(path, files(trc_file_index));
        experiment = experiment{1};
        has_2kg_weight = str2num(experiment(end-5));      % based on file name
    else
        experiment = append(path,files);
        has_2kg_weight = str2num(experiment(end-5));      % based on file name
    end
    
    % consider the correct model in the analysis, based on the .trc files
    
    [aux_optimization_status, aux_unfeasibility_flags, tOptim(trc_file_index), aux_result_file] = RMR_analysis(dataset_considered, model_2kg, ...
        experiment, 0, weight_coord, time_interval, dynamic_bounds, ...
        enforce_GH_constraint, execlude_locked,execl_clav, execl_wing, ...
        t_end, saving_path, tag_con, tag_cost);

    optimizationStatus(trc_file_index).experiment = aux_optimization_status;
    result_file_RMR{trc_file_index} = aux_result_file;
    unfeasibility_flag(trc_file_index).experiment = aux_unfeasibility_flags;
    fprintf('\n Solved with %i unfeasible solutions \n \n \n', sum(aux_unfeasibility_flags));
end

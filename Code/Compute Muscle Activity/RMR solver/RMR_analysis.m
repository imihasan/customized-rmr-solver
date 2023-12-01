function [optimizationStatus, unfeasibility_flags, tOptim, file_results] = RMR_analysis(subject_considered, ...
    model_original, trc_file, motion_file, weight_coord, time_interval, dynamic_activation_bounds, ...
    flag_GH_enforced, flag_constraint, flag_excludeSpine, flag_excludeNeck, flag_excludeBoth,execlude_locked,saving_path)
% Rapid Muscle Redundancy (RMR) solver, leveraging OpenSim API.
% Starting from experimental marker data (in .trc format) the optimal
% muscle activations are found that can reproduce the motion, solving:
%
%   min   sum (w_i * a_i^2)+ sum (w_j * c_j^2)
%   a,c    i                  j                 
%
%   s.t.        a_min<=a_i<=a_max    for muscle activations
%              -l<=c_j<=l            for coordinateActuators controls (if present)
%        acc_{j,FD} = acc_{j,data}   constraint on accelerations
%               F_{GH} \in Cone      glenohumeral constraint
%
%
% The code is written specifically to consider a thoracoscapular shoulder
% model that has been already scaled to the biometrics of the subject of
% interest However, this script can be generalized to consider other models 
% and data without changing its main structure. 
%
% INPUTS:
% * subject_considered: string defining the name of the subject analyzed
%                       (used to save results)
% * model_original: model to be used for the analysis
%                  (valid TSM model with GH markers and coordinate actuators)
% * trc_file : path to the file - and file name - from which to retrieve
%              marker data for IK and subsequent RMR analysis (set to 0 if the
%              input is actually the motion file)
% * motion_file: path to the file - and file name - that carries
%                information on the coordinates (set to 0 if trc file is
%                used)
% * weight_coord : 4x1 vector indicating the weight of each scapula DoF in
%                  the IK tracking 
%                  (order: abduction, elevation, upward rotation, winging)
% * time_interval : downsampling of the original data, to reduce
%                   computation effort for the RMR. For example, if set to 10, 
%                   every 10th time point is selected.
% *dynamic_activation_bounds : flag to indicate whether dynamic bounds must
%                              be used to limit the activation values during 
%                              RMR solution
% * flag_GH_enforced: true or false, if glenohumeral constraint is
%                     considered or not
% * saving_path: path to where the results of the redundancy solver are
%                saved
%
% OUTPUT:
% * optimizationStatus: struct containing the status of the optimization at
%                       each of the timestep in which RMR was performed.
% * unfeasibility_flags: an array of the same length as the time instants
%                        considered, containing 0 if the problem is solved and  
%                        1 if it is unfeasible
% * tOptim : time required to perform the complete optimization
% * file_results: path and name of the file where the activation results
%                 are saved
% The function also saves plots of the analysis performed, and the muscle
% activation variables (together with coordinate actuators controls) in a
% .mat file
% author: Italo Belli (i.belli@tudelft.nl) 2022
%Edited By: Ibrahim Mohammed Hasan (imihasan@kth.se) 2023, KTH MoveAbility Lab,
%KTH Royal Institute of Technology, Stockholm, Sweden.

%% Import the OpenSim libraries.
import org.opensim.modeling.*;

%% General settings
% if these are set to true, results are printed but the code will be slower
print_flag = true;         
withviz = false;

%% Set the correct paths
% set the path current folder to be the one where this script is contained
mfile_name          = mfilename('fullpath');
[pathstr,~,~]  = fileparts(mfile_name);
cd(pathstr);

% getting path to other folders in this repo
addpath(pathstr)
cd ..\..\..\
path_to_repo = pwd;
addpath(path_to_repo)
addpath(fullfile(path_to_repo, 'Code\Data Processing\'))

% cd to Personal Results to have all the results saved there
cd(saving_path);

% create a temporary copy of the model, to be used in the function. In this
% way, the model can be modified freely here without interfering with its
% state/properties outside this function
model_temp = model_original.clone();     

%% Getting quantities about GlenoHumeral joint
% get the glenohumeral joint
alljoints = model_temp.getJointSet;
glen = alljoints.get('GlenoHumeral');
Lglen = alljoints.get('LGlenoHumeral'); %Read the left Glenohumeral Joint As well

state = model_temp.initSystem();
[maxAngle, LmaxAngle,~,~] = get_glenoid_status(model_temp, state); % the value for maxAngle can also be given directly by the user

%% Load the trc file to be considered, if the input is a trc file, and perform IK
if trc_file
    [~, experiment_name] = fileparts(trc_file);
    [markersExp, timesExp, ~, unitsExp] = readTRC(trc_file);
    start_time = timesExp(1);
    end_time =  timesExp(end);
    
    if strcmp(unitsExp, 'mm')
        markersExp = markersExp/1000;
        unitsExp = 'm';
    end
    
    frequency_trc_data = 1/(timesExp(2)-timesExp(1));

    % getting the values of default scapula coordinate 
    % we get the values of the coordinates describing the scapula position from 
    % the general model in default pose
    %Right
    scapula_abd = model_temp.getJointSet().get("scapulothoracic").get_coordinates(0);
    scapula_ele = model_temp.getJointSet().get("scapulothoracic").get_coordinates(1);
    scapula_urt = model_temp.getJointSet().get("scapulothoracic").get_coordinates(2);
    scapula_wng = model_temp.getJointSet().get("scapulothoracic").get_coordinates(3);

    default_sa = scapula_abd.get_default_value();
    default_se = scapula_ele.get_default_value();
    default_su = scapula_urt.get_default_value();
    default_sw = scapula_wng.get_default_value();

    %Left
    Lscapula_abd = model_temp.getJointSet().get("Lscapulothoracic").get_coordinates(0);
    Lscapula_ele = model_temp.getJointSet().get("Lscapulothoracic").get_coordinates(1);
    Lscapula_urt = model_temp.getJointSet().get("Lscapulothoracic").get_coordinates(2);
    Lscapula_wng = model_temp.getJointSet().get("Lscapulothoracic").get_coordinates(3);

    Ldefault_sa = Lscapula_abd.get_default_value();
    Ldefault_se = Lscapula_ele.get_default_value();
    Ldefault_su = Lscapula_urt.get_default_value();
    Ldefault_sw = Lscapula_wng.get_default_value();

    
    % Performing IK
    % perform IK on the basis of marker data to retrieve the motion file for
    % the coordinates of the model
    
    motion_file_name = append(experiment_name, '.mot');
    
    ikSetupFile = fullfile(path_to_repo,'ExperimentalData', 'IK setup files', 'IKSetup.xml');
    
    ikTool = InverseKinematicsTool(ikSetupFile);
    ikTool.setMarkerDataFileName(trc_file);
    ikTool.setOutputMotionFileName(fullfile(saving_path, motion_file_name));
    ikTool.set_report_marker_locations(1);
    ikTool.setStartTime(start_time);
    %The end time is set to quarter of the total time. That's to run a
    %short inverse kinmatics to get the initial values of the ground pelvis
    %coordinates to position the model correctly by change the default
    %values of these coordinates
    ikTool.setEndTime(end_time/4);                                          
    ikTool.setModel(model_temp);
    
    % set the reference values for the scapula coordinates (last 4 tasks)
    %num_IK_tasks = ikTool.getIKTaskSet.getSize(); %No need for that, tasks are called by their names instead
    
    %set the weight of each coordinate in the tracking tasks
    % Right
    ikTool.getIKTaskSet.get("scapula_abduction").setWeight(weight_coord(1));
    ikTool.getIKTaskSet.get("scapula_elevation").setWeight(weight_coord(2));
    ikTool.getIKTaskSet.get("scapula_upward_rot").setWeight(weight_coord(3));
    ikTool.getIKTaskSet.get("scapula_winging").setWeight(weight_coord(4));

    %Left
    ikTool.getIKTaskSet.get("Lscapula_abduction").setWeight(weight_coord(1));
    ikTool.getIKTaskSet.get("Lscapula_elevation").setWeight(weight_coord(2));
    ikTool.getIKTaskSet.get("Lscapula_upward_rot").setWeight(weight_coord(3));
    ikTool.getIKTaskSet.get("Lscapula_winging").setWeight(weight_coord(4));

    % set also the values here
    %Right
    IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get("scapula_abduction")).setValue(default_sa);
    IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get("scapula_elevation")).setValue(default_se);
    IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get("scapula_upward_rot")).setValue(default_su);
    IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get("scapula_winging")).setValue(default_sw);
    % 
    %Left
    IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get("Lscapula_abduction")).setValue(Ldefault_sa);
    IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get("Lscapula_elevation")).setValue(Ldefault_se);
    IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get("Lscapula_upward_rot")).setValue(Ldefault_su);
    IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet.get("Lscapula_winging")).setValue(Ldefault_sw);
    ikTool.print(fullfile(path_to_repo,'ExperimentalData', 'IK setup files','RMR_autogenerated_IK_setup.xml'));
    
    ikTool.run();

    %Get Coordinates to adjust initial position of the model. The motion
    %file is imported and the mean value of the coordinates of interest is
    %obtained 
    lowpassFreq = 3.0; % Hz
    timeRange = [start_time end_time];
    [coordinates, coordNames, ~] = loadFilterCropArray(motion_file_name, lowpassFreq, timeRange); 
    
    ground_coord=["ground_pelvis_coord_3", "ground_pelvis_coord_4",  "ground_pelvis_coord_5"]; % Coordinates that will be adjusted

    for gcord=1:length( ground_coord)
        defval=mean(coordinates(:,find(coordNames==ground_coord(gcord))));  %Get the mean value of such coordinates
        model_temp.updCoordinateSet().get(ground_coord(gcord)).setDefaultValue(defval); %Set the default value
        model_temp.updCoordinateSet().get(ground_coord(gcord)).set_locked(1); %Lock the coordinates on interest
    end
    % model_temp.finalizeConnections();
    % model_temp.print(fullfile(saving_path,"model.osim"));
    %Run Inverse Kinematics again but for the full time range
    ikTool.setModel(model_temp); %Set the model template after adjusting the coordinates
    ikTool.setEndTime(end_time); %Set the end time to cover full range
    ikTool.run();

else 
    [~, experiment_name] = fileparts(motion_file);
    motion_file_name = motion_file;
    q = read_motionFile(motion_file_name);
    time = q.data(:,1); 
    start_time = time(1);
    end_time = time(end);
    frequency_trc_data = 1/(time(2)-time(1));
end


%% getting the kinematic data that we need
% Use the loadFilterCropArray() function provided by OpenSim Tutorial to load the 
% coordinate kinematic and generalized force data into MATLAB arrays. This 
% function also filters and crops the loaded array based on its two input 
% arguments (more details in loadFilterCropArray.m).
% model_temp=Model(fullfile(saving_path,"model.osim"));



lowpassFreq = 3.0; % Hz
timeRange = [start_time end_time];

% get the coordinates from the output of the IK in rad for the rotational
% joints
[coordinates, coordNames, timesExp] = loadFilterCropArray(motion_file_name, lowpassFreq, timeRange); 


%Drop Locked Coordinates
if execlude_locked==1
    lc=1;
    for acc=1:size(coordNames)
        if (model_temp.getCoordinateSet().get(coordNames(acc)).get_locked())==1
            indx_coord(lc)=acc;
            locked_acts_names(lc)=string(model_temp.getCoordinateSet().get(coordNames(acc)).getName())+"_actuator";
            lc=lc+1;
        end
    end
    %delete locked coordinates
    coordinates(:,indx_coord)=[];
    coordNames(indx_coord)=[];
end


%Exclude coordinates that are not necessary in the solving the muscles
%redundancy
if flag_excludeSpine==1
    ex_Names=["L5_S1_Flex_Ext", "L1_T12_axial_rotation"];
    indx1=find(coordNames==ex_Names(1));
    indx2=find(coordNames==ex_Names(2));
    coordinates(:,indx1:indx2)=[];
    coordNames(indx1:indx2)=[];
elseif flag_excludeNeck==1
    ex_Names=["T1_C7_Flex_Ext", "C2_C1_axial_rotation"];
    indx1=find(coordNames==ex_Names(1));
    indx2=find(coordNames==ex_Names(2));
    coordinates(:,indx1:indx2)=[];
    coordNames(indx1:indx2)=[];
elseif flag_excludeBoth==1
    %Drop Spine Coordinates
    ex_Names_S=["L5_S1_Flex_Ext", "L1_T12_axial_rotation"];
    indx1=find(coordNames==ex_Names_S(1));
    indx2=find(coordNames==ex_Names_S(2));
    coordinates(:,indx1:indx2)=[];
    coordNames(indx1:indx2)=[];

    %Drop Neck Coordinates
    ex_Names_N=["T1_C7_Flex_Ext", "C2_C1_axial_rotation"];
    indx1=find(coordNames==ex_Names_N(1));
    indx2=find(coordNames==ex_Names_N(2));
    coordinates(:,indx1:indx2)=[];
    coordNames(indx1:indx2)=[];
end


%Convert rotational degrees from degrees to radians
trans_coord=["GH_Tx", "GH_Ty", "GH_Tz", "LGH_Tx", "LGH_Ty", "LGH_Tz", "ground_pelvis_coord_3", "ground_pelvis_coord_4",  "ground_pelvis_coord_5"];

% coordinates=deg2rad(coordinates);
for coor=1:size(coordNames)
    if ~ismember(coordNames(coor),trans_coord)
        coordinates(:,coor)=deg2rad(coordinates(:,coor));
    end
end

% get the velocities for each joint in rad/s
time_step_data = timesExp(2)-timesExp(1);
speeds = zeros(size(coordinates));
for i=1:size(coordNames,1)
    speeds(:,i) = gradient(coordinates(:,i), time_step_data);
end
speedNames = coordNames;

% get the accelerations for each coordinate in rad/s^2
accelerations = zeros(size(speeds));
for i=1:size(coordNames,1)
    accelerations(:,i) = gradient(speeds(:,i), time_step_data);
end
accNames = speedNames;




% visually check the values of joint states, speeds and accelerations
if print_flag
    Gr = ceil(sqrt(length(coordNames)));
    figure
    for i=1:length(coordNames)
    subplot(Gr,Gr,i)
    hold on
    plot(coordinates(:,i))
    plot(speeds(:,i))
    plot(accelerations(:,i))
    title(coordNames{i});
    hold off
    grid on
    end
    legend("coords", "speeds", "accs")
end

%% Store max isometric force values and disable muscle dynamics
muscles = model_temp.getMuscles();
numMuscles = muscles.getSize();
muscles_downcasted = cell(numMuscles,1);
muscleNames = cell(numMuscles,1);

% save here downcasted muscles to a list 
for index_muscle = 1:numMuscles
   % Downcast base muscle to Millard2012EquilibriumMuscle
   muscles_downcasted(index_muscle) = Millard2012EquilibriumMuscle.safeDownCast(muscles.get(index_muscle-1));
   muscleNames{index_muscle} = char(muscles_downcasted{index_muscle});
   muscles_downcasted{index_muscle}.set_ignore_tendon_compliance(true);  % not really relevant as actuation will be overwritten
   muscles_downcasted{index_muscle}.set_ignore_activation_dynamics(true);
end

if (withviz == true)
    model_temp.setUseVisualizer(true);
end

% Update the system to include any muscle modeling changes
state = model_temp.initSystem();

%% Get coordinate actuators 
allActs = model_temp.getActuators; 
num_acts = getSize(allActs); 
acts = cell(num_acts,1);
actsNames=cell(num_acts,1);
% get all actuators and override actuation for the muscles only
for i = 1:num_acts
    acts(i) = ScalarActuator.safeDownCast(allActs.get(i-1));
    actsNames(i)=allActs.get(i-1).getName();
    if i<=numMuscles
        acts{i}.overrideActuation(state, true); %Enable Muscles Actuation (i.e. apply force)
    end
end

%Exclude locked coordinate actuators
if execlude_locked
    unloc=1;
    for unloc_act=1:num_acts
        if ~ismember(actsNames(unloc_act),locked_acts_names)
            acts_unlock(unloc)=acts(unloc_act);
            actsNames_unloc(unloc)=actsNames(unloc_act);
            unloc=unloc+1;
        end
    end
    acts=acts_unlock;
    actsNames=actsNames_unloc;
end

%exclude coordinate actuators that are not necessary in the solving the muscles
%redundancy
if flag_excludeSpine==1
    ex_Names=["L5_S1_Flex_Ext_actuator", "L1_T12_axial_rotation_actuator"];
    indx1=find(actsNames==ex_Names(1));
    indx2=find(actsNames==ex_Names(2));
    acts(indx1:indx2)=[];
    actsNames(indx1:indx2)=[];
elseif flag_excludeNeck==1
    ex_Names=["T1_C7_Flex_Ext_actuator", "C2_C1_axial_rotation_actuator"];
    indx1=find(actsNames==ex_Names(1));
    indx2=find(actsNames==ex_Names(2));
    acts(indx1:indx2)=[];
    actsNames(indx1:indx2)=[];
elseif flag_excludeBoth==1

    ex_Names_S=["L5_S1_Flex_Ext_actuator", "L1_T12_axial_rotation_actuator"];
    indx1=find(actsNames==ex_Names_S(1));
    indx2=find(actsNames==ex_Names_S(2));
    acts(indx1:indx2)=[];
    actsNames(indx1:indx2)=[];

    ex_Names_N=["T1_C7_Flex_Ext_actuator", "C2_C1_axial_rotation_actuator"];
    indx1=find(actsNames==ex_Names_N(1));
    indx2=find(actsNames==ex_Names_N(2));
    acts(indx1:indx2)=[];
    actsNames(indx1:indx2)=[];
end

%Recall the new size of the actuators
num_acts=length(actsNames);
%% Deactivate Constraints 
%Constrained are deactivated in the dynamic analysis to avoid generalized
%forces imposed by which are not physical
if flag_constraint==1
    for cons=1:model_temp.updConstraintSet().getSize()
        model_temp.updConstraintSet().get(cons-1).set_isEnforced(0);
    end
end
%% Perform optimization
% We use FMINCON to solve the optimization problem at selected time points. 
% The 'time_interval' variable selects the time points to be included in the 
% optimization. For example, if set to 10, every 10th time point is selected. A 
% time interval of 1 will select all available time points. 
time_step_RMR = time_step_data * time_interval;

% Update data arrays based on the time_interval.
N = size(coordinates, 1);
coordinates = coordinates(1:time_interval:N, :);
speeds = speeds(1:time_interval:N, :);
accelerations = accelerations(1:time_interval:N, :);
numTimePoints = size(coordinates, 1);
unfeasibility_flags = zeros(size(numTimePoints));

% Create the FMINCON options structure.
options = optimoptions('fmincon','Display','none', ...
     'TolCon',1e-2,'TolFun',1e-3,'TolX',1e-3,'MaxFunEvals',100000, ...
     'MaxIter',10000,'Algorithm','sqp', 'StepTolerance', 1e-8); %, 'DiffMinChange', 1.0e-2);
 
% Construct initial guess and bounds arrays
numCoords = length(coordNames);
numCoordActs = num_acts-numMuscles;
k = inf;
t_act = 0.01;           % activation time constant for muscles
t_deact = 0.04;         % deactivation time constant

lb = [zeros(1,numMuscles), -k*ones(1,numCoordActs)];
ub = [ones(1,numMuscles), k*ones(1,numCoordActs)];
x0 = [0.1* ones(1,numMuscles), zeros(1,numCoordActs)];

% We define the activation squared cost as a MATLAB anonymous function
% It is model specific!
epsilon = 0;
%The Weighting Matrix needs to change. Dimensions are not consistent
%w = [ones(1,numMuscles), epsilon*ones(1,8), 10*ones(1,9)];     % the cost function is written such that it allows the use of coord acts for the underactuated coordinates
%Now the weights are consistent


joints_penalized=["elbow","GlenoHumeral","scapulothoracic", "radioulnar",
    "Lelbow","LGlenoHumeral","Lscapulothoracic", "Lradioulnar"];
wc=zeros(1,numCoordActs);
for w_coord=1:length(coordNames)
    jnt_name=string(model_temp.getCoordinateSet().get(coordNames(w_coord)).getJoint().getName());
    if ismember(jnt_name,joints_penalized)
        wc(w_coord)=10;
    else
        wc(w_coord)=0;
    end
end

w = [ones(1,numMuscles), wc];
cost =@(x) sum(w.*(x.^2));

% Pre-allocate arrays to be filled in the optimization loop
fl = zeros(1, numMuscles);
fv = zeros(1, numMuscles);
fp = zeros(1, numMuscles);
cosPenn = zeros(1, numMuscles);
Fmax = zeros(1, numMuscles);
A_eq_acc = zeros(numCoords,num_acts);
% LA_eq_acc = zeros(numCoords,num_acts); %Left

A_eq_force = zeros(3, num_acts);
LA_eq_force = zeros(3, num_acts); %Left

xsol = zeros(numTimePoints, length(x0));
simulatedAccelerations = zeros(numTimePoints, length(coordNames));
optimizationStatus = cell(numTimePoints,1);

norm_fv_in_ground = zeros(numTimePoints, 3);
Lnorm_fv_in_ground = zeros(numTimePoints, 3); %Left

norm_fv_rotated = zeros(numTimePoints, 3);
Lnorm_fv_rotated = zeros(numTimePoints, 3); %Left

rel_angle = zeros(numTimePoints,1);
Lrel_angle = zeros(numTimePoints,1); %Left

%Initialize variables for storing ligament forces
ligaments=["s_glenohum_r"; "m_glenohum_r";  "i_glenohum_r"; "coracohum_r";
    "s_glenohum_l"; "m_glenohum_l";  "i_glenohum_l"; "coracohum_l"]; %Ligaments included in the model
lig_force=zeros(numTimePoints,length(ligaments)); %initialize the force vector where ligaments forces will be stored
lig_ratio=zeros(numTimePoints,length(ligaments));

% get model quantities we still need
coords = model_temp.getCoordinateSet();

for index_muscle = 1:numMuscles
    Fmax(index_muscle) = muscles_downcasted{index_muscle}.getMaxIsometricForce();
end

% do not track plane_elv and axial_rot during shrugging, as they are poorly 
% defined when humerus is vertical. Similar in what done 
% by Seth et al. in https://simtk.org/projects/thoracoscapular we lock these
if strcmpi(experiment_name(1:5), 'shrug')
    model_temp.getCoordinateSet().get('plane_elv').set_default_value(-0.433725);
    model_temp.getCoordinateSet().get('plane_elv').set_locked(true);
    model_temp.getCoordinateSet().get('axial_rot').set_default_value(0.8125346);
    model_temp.getCoordinateSet().get('axial_rot').set_locked(true);
end

tb_Optim=zeros(1,length(timesExp));
tA_Optim=zeros(1,length(timesExp));
tf_Optim=zeros(1,length(timesExp));

% create a struct containing relevant information to be passed to the
% function simulating the accelerations and reaction forces and moments
% induced in the model
%Those parameters were defined inside the loop and were taking
%computational time specially params.Lglen = Lglen I moved it outside to
%save compuational time and kept only parmeters that needs to be updted
%inside the loop
params.model = model_temp; 
params.state = state;    
params.coords = coords; 
params.coordNames = coordNames; 
params.acts = acts; 
params.muscles = muscles; 
params.numMuscles = numMuscles; 
params.useMuscles = 1; 
params.useControls = 1; 
params.glen = glen; 
params.Lglen = Lglen; 

tic
% enter in the optimization loop
for time_instant = 1:numTimePoints
    tb_tic=tic;
    if print_flag
        fprintf('.');
        if(mod(time_instant,80) == 0)
            fprintf('\n %i', time_instant);
        end
    end

    % set the time of the simulation to be the current one
    % (this is especially important if an external force is present, so
    % that the force value is applied at the right instant of time)
    state.setTime((time_instant-1)*time_interval/frequency_trc_data)

    % Loop through model coordinates to set coordinate values and speeds. We set
    % all coordinates to make sure we have the correct kinematic state when 
    % compute muscle multipliers and moment arms.
    for j = 1:length(coordNames)
        coord = coords.get(coordNames{j});
        coord.setValue(state, coordinates(time_instant,j), false); % instead of fals replace so that does the assembly on the last call (j==length)
        coord.setSpeedValue(state, speeds(time_instant,j));
    end

    % realize the system to the velocity stage
    model_temp.realizeVelocity(state);
   
    
    % equilibrate the muscles to make them start in the correct state
    model_temp.equilibrateMuscles(state);
    
    modelControls = model_temp.getControls(state);

    model_temp.realizeAcceleration(state); 
    for g=1:length(ligaments)
        lig=Ligament.safeDownCast(model_temp.updForceSet().get(ligaments(g)));
        lig_ratio(time_instant,g)=lig.getLength(state)/lig.getRestingLength(); %Tension=0 if ratio<=0
        lig_force(time_instant,g)=lig.getTension(state);
    end

    % Populate the muscle multiplier arrays. To do this, we must have realized 
    % the system to the velocity stage
    for index_muscle = 1:numMuscles
        fl(index_muscle) = muscles_downcasted{index_muscle}.getActiveForceLengthMultiplier(state);            
        fv(index_muscle) = muscles_downcasted{index_muscle}.getForceVelocityMultiplier(state);
        fp(index_muscle) = muscles_downcasted{index_muscle}.getPassiveForceMultiplier(state);
        cosPenn(index_muscle) = muscles_downcasted{index_muscle}.getCosPennationAngle(state);
    end 

    % get the vector Vec_H2GC between humeral head and the glenoid center
    % (it is expressed in the ground frame)
    [~,~,Vec_H2GC, LVec_H2GC] = get_glenoid_status(model_temp, state); %Left LVec_H2GC is called as output
    
    % store the values of active and passive maximum force in the current
    % configuration
    AMuscForce = (fl.*fv.*Fmax.*cosPenn)'; 
    PMuscForce = (Fmax.*fp.*cosPenn)'; 
    

    %Define the remaining params structs elements
    params.AMuscForce = AMuscForce;
    params.PMuscForce = PMuscForce;
    params.modelControls = modelControls;


   
    [q_ddot_0, F_r0, LF_r0,~,~] = findInducedAccelerationsForceMomentsGH(zeros(1,num_acts), params); %Get Left Forces in the output
    delQ_delX = eye(num_acts);

    for k = 1:num_acts
        [incrementalForceAccel_k, F_rk, LF_rk, ~, ~] = findInducedAccelerationsForceMomentsGH(delQ_delX(k,:),params); %Get Left Forces in the output
        kthColumn_A_eq_acc =  incrementalForceAccel_k - q_ddot_0;
        A_eq_acc(:,k) = kthColumn_A_eq_acc;

        %Include the left side here too
        kthColumn_A_eq_force =  F_rk - F_r0;
        LkthColumn_A_eq_force =  LF_rk - LF_r0; %Left Side
        A_eq_force(:,k) = kthColumn_A_eq_force;
        LA_eq_force(:,k) = LkthColumn_A_eq_force; %Left Side
    end

    Beq = accelerations(time_instant,:)' - q_ddot_0;

    % if the task considered is a SHRUGGING task, do not track the 'plane
    % of elevation' (13th) and the 'axial rotation' (15th) coordinates as 
    % they are poorly defined.
    % if strcmpi(experiment_name(1:5), 'shrug')
    %     A_eq_acc(15, :) = zeros(size(A_eq_acc(15, :)));
    %     A_eq_acc(13, :) = zeros(size(A_eq_acc(13, :)));
    %     Beq(15, :) = zeros(size(Beq(15, :)));
    %     Beq(13, :) = zeros(size(Beq(13, :)));
    % end

    %Get Ligaments forces
    %This part has been added to realize ligament forces at each time step
    %and export this in the solution

    tb_Optim(time_instant)=toc(tb_tic); %Measure time taken until entering the fmincon function

    tA_tic=tic;
    % Call FMINCON to solve the problem
    if flag_GH_enforced
        [x,~,exitflag,output] = fmincon(cost, x0, [], [], A_eq_acc, Beq, lb, ub, @(x)jntrxncon_linForce(x, Vec_H2GC, LVec_H2GC, maxAngle, LmaxAngle, A_eq_force, LA_eq_force, F_r0, LF_r0,1), options);
        if exitflag ==0
            % call the solver again, starting from current x, in case the maximum iterations are exceeded
            [x,~,exitflag,output] = fmincon(cost, x, [], [], A_eq_acc, Beq, lb, ub, @(x)jntrxncon_linForce(x, Vec_H2GC, LVec_H2GC, maxAngle, LmaxAngle, A_eq_force, LA_eq_force, F_r0, LF_r0,1), options);
        end
        if exitflag<0 && time_instant>1
            % call the solver again, starting from previous optimum found,
            % in case optimization gets stuck in local minimum 
            [x,~,exitflag,output] = fmincon(cost, xsol(time_instant-1, :), [], [], A_eq_acc, Beq, lb, ub, @(x)jntrxncon_linForce(x, Vec_H2GC, LVec_H2GC, maxAngle, LmaxAngle, A_eq_force, LA_eq_force, F_r0, LF_r0,1), options);
        end
    else
        [x,~,exitflag,output] = fmincon(cost, x0, [], [], A_eq_acc, Beq, lb, ub, [], options);
        if exitflag ==0
            % call the solver again, starting from current x, in case the maximum iterations are exceeded
            [x,~,exitflag,output] = fmincon(cost, x, [], [], A_eq_acc, Beq, lb, ub, [], options);
        end
        if exitflag<0 && time_instant>1
            % call the solver again, starting from previous optimum found,
            % in case optimization gets stuck in local minimum 
            [x,~,exitflag,output] = fmincon(cost, xsol(time_instant-1, :), [], [], A_eq_acc, Beq, lb, ub, [], options);
        end
    end
    tA_Optim(time_instant)=toc(tA_tic);
    optimizationStatus{time_instant} = output;

    tf_tic=tic;
    if exitflag<0
        unfeasibility_flags(time_instant) = 1;
    end

    % get best feasible point, if different from what returned by fmincon
    if size(output.bestfeasible,1)>0
        x = output.bestfeasible.x;
    end
    
    % Store solution
    xsol(time_instant, :) = x;

    % dynamically update the upper and lower bounds for the activations
    if dynamic_activation_bounds
        for k = 1:numMuscles
            lb(k) = max(x(k) - x(k) * (0.5 + 1.5 * x(k)) * time_step_RMR /t_deact, 0);
            ub(k) = min (x(k) + (1-x(k)) * time_step_RMR / (t_act * (0.5 + 1.5*x(k))), 1);
        end
    end

    % if we want to print suff, we need to compute it now
    if print_flag
        % Retrieve the optimal accelerations
        [simulatedAccelerations(time_instant,:), ~, ~] = findInducedAccelerationsForceMomentsGH(x,params);

        % if the task considered is a SHRUGGING task, disregard
        % accelerations for plane of elevation and axial rotation (that are
        % locked)
        if strcmpi(experiment_name(1:5), 'shrug')
            simulatedAccelerations(time_instant,[13,15]) = 0;
        end

        % retrieve the position of the joint reaction force on the approximated
        % glenoid computing the reaction force vector at the given joint
        % The force is expressed in the ground frame
        %Right
        force_vec = A_eq_force * xsol(time_instant, :)' + F_r0;
        force_vecc(time_instant,:)=A_eq_force * xsol(time_instant, :)' + F_r0;
        %Left
        Lforce_vec = LA_eq_force * xsol(time_instant, :)' + LF_r0;
        Lforce_vecc(time_instant,:) = LA_eq_force * xsol(time_instant, :)' + LF_r0;
        
        % evaluate the relative angle between the reaction force and Vec_H2GC
        %Right
        cosTheta = max(min(dot(Vec_H2GC,force_vec)/(norm(Vec_H2GC)*norm(force_vec)),1),-1);
        rel_angle(time_instant) = real(acosd(cosTheta));
        %Left
        LcosTheta = max(min(dot(LVec_H2GC,Lforce_vec)/(norm(LVec_H2GC)*norm(Lforce_vec)),1),-1);
        Lrel_angle(time_instant) = real(acosd(LcosTheta));
    
        % evaluate the position on the glenoid where reaction force is exerted
        %Right
        norm_Vec_H2GC = Vec_H2GC/norm(Vec_H2GC);
        norm_fv_in_ground(time_instant,:) = force_vec/norm(force_vec);
        %Left
        Lnorm_Vec_H2GC = LVec_H2GC/norm(LVec_H2GC);
        Lnorm_fv_in_ground(time_instant,:) = Lforce_vec/norm(Lforce_vec);

        %Right
        beta_angle = atan(norm_Vec_H2GC(3)/norm_Vec_H2GC(1));
        alpha_angle = atan(norm_Vec_H2GC(3)/(sin(beta_angle)*norm_Vec_H2GC(2)));
        %Left
        Lbeta_angle = atan(Lnorm_Vec_H2GC(3)/Lnorm_Vec_H2GC(1));
        Lalpha_angle = atan(Lnorm_Vec_H2GC(3)/(sin(Lbeta_angle)*Lnorm_Vec_H2GC(2)));

        %Intrinsic rotation sequence around y and z respectively to align the global
        %frame with the glenoid cavity frame whose z axis is -norm_Vec_H2GC
        %Right
        Ry = [cos(beta_angle) 0 sin(beta_angle); 0 1 0; -sin(beta_angle) 0 cos(beta_angle)];
        Rz = [cos(alpha_angle) -sin(alpha_angle) 0; sin(alpha_angle) cos(alpha_angle) 0; 0 0 1];
        %Left
        LRy = [cos(Lbeta_angle) 0 sin(Lbeta_angle); 0 1 0; -sin(Lbeta_angle) 0 cos(Lbeta_angle)];
        LRz = [cos(Lalpha_angle) -sin(Lalpha_angle) 0; sin(Lalpha_angle) cos(Lalpha_angle) 0; 0 0 1];
    
        
        %Right 
        norm_fv_rotated(time_instant,:) = Rz*Ry*norm_fv_in_ground(time_instant,:)';
        force_vecc_rotated(time_instant,:)=Rz*Ry*force_vecc(time_instant,:)';
        %Left 
        Lnorm_fv_rotated(time_instant,:) = LRz*LRy*Lnorm_fv_in_ground(time_instant,:)';
        Lforce_vecc_rotated(time_instant,:)=Rz*Ry*Lforce_vecc(time_instant,:)';
    end

    if (withviz == true)
        model_temp.getVisualizer.show(state);
    end
    tf_Optim(time_instant)=toc(tf_tic);
end

tOptim = toc;

%% Plot results
% According to the value of the 'print_flag'
if print_flag
    % plot muscle activations
    f1 = figure;
    title("Muscle Activations")
    muscleNames = ArrayStr();
    muscles.getNames(muscleNames);
    Dim = ceil(sqrt(numMuscles)); %Ajusted the the subplot Array
    pgc = linspace(0, 100, size(xsol,1));
    for i = 1:numMuscles
       subplot(Dim,Dim,i)
       hold on
       plot(pgc,xsol(:,i),'b-')
       ylim([0 1])
       muscName = muscleNames.get(i-1).toCharArray';
       title(muscName(1:end), 'interpreter', 'none')
       hold off
    end
    legend("muscle activation")
    f1.WindowState = 'maximized';
    name_fig1 = append(experiment_name, '_MuscleActivations.png');
    name_fig1svg = append(experiment_name, '_MuscleActivations.svg');
    saveas(f1, name_fig1)
    saveas(f1, name_fig1svg)
    
    % Plot reserve actuator excitations.
    f2 = figure;
    title("Reserve actuators")
    side = ceil(sqrt(numCoordActs));
    for i = 1:numCoordActs
        subplot(side,side,i)
        hold on
        plot(pgc, xsol(:,numMuscles+i), 'linewidth', 2);
        title(char(acts{numMuscles+i}),'Interpreter','none');
        hold off
    end
    legend("reserve act value")
    f2.WindowState = 'maximized';
    name_fig2 = append(experiment_name, '_ReserveActuators.png');
    name_fig2svg = append(experiment_name, '_ReserveActuators.svg');
    saveas(f2, name_fig2)
    saveas(f2, name_fig2svg)

    % plot accelerations
    f3 = figure;
    title("Accelerations")
    side = ceil(sqrt(length(coordNames)));
    for i = 1:length(coordNames)
        subplot(side,side,i)
        hold on
        plot(accelerations(:, i), 'linewidth', 1.5);
        plot(simulatedAccelerations(:, i), 'linewidth', 1.5);
        xlabel("samples")
        ylabel("[]/s^2")
        grid on
        title(coordNames{i},'Interpreter','none');
        hold off
    end
    legend("measured", "simulated")
    f3.WindowState = 'maximized';
    name_fig3 = append(experiment_name, '_AccelerationsMatching.png');
    name_fig3svg = append(experiment_name, '_AccelerationsMatching.svg');
    saveas(f3, name_fig3)
    saveas(f3, name_fig3svg)

    % plot the constraint violation on the accelerations per coordinate
    violation = abs(accelerations-simulatedAccelerations);

     if strcmpi(experiment_name(1:5), 'shrug')
        violation(:,[13,15]) = 0;
     end

    f4 = figure;
    for i = 1:length(coordNames)
        subplot(side,side,i)
        hold on
        plot(violation(:,i), 'linewidth', 1.5);
        xlabel("samples")
        ylabel("[]/s^2")
        grid on
        title(coordNames{i},'Interpreter','none');
        hold off
    end
    legend("acc violation")
    f4.WindowState = 'maximized';
    name_fig4 = append(experiment_name, '_AccViolation.png');
    name_fig4svg = append(experiment_name, '_AccViolation.svg');
    saveas(f4, name_fig4)
    saveas(f4, name_fig4svg)

    % plot the constraint violation per timestep
    violation_t = sum(violation, 2);
    f5 = figure;
    hold on
    scatter(1:numTimePoints ,violation_t, 'filled')
    plot(1:numTimePoints, violation_t, 'blue')
    xlabel("samples")
    ylabel("const violation")
    grid on
    title("Cumulative constraint violation per time-step")
    hold off
    f5.WindowState = 'maximized';
    name_fig5 = append(experiment_name, '_CumulativeAccViolation.png');
    name_fig5svg = append(experiment_name, '_CumulativeAccViolation.svg');
    saveas(f5, name_fig5)
    saveas(f5, name_fig5svg)

    % plot the position of the GH force on the glenoid
    radius = sind(maxAngle);
    p=nsidedpoly(1000, 'Center', [0,0], 'Radius', radius);
    c = linspace(0,timesExp(end),length(norm_fv_rotated));
    f6 = figure;
    hold on
    plot(p, 'FaceColor', 'r')
    for time_instant=1:numTimePoints
        scatter(-norm_fv_rotated(time_instant,3), -norm_fv_rotated(time_instant,1), [], c(time_instant), 'filled')
    end
    hcb = colorbar;
    h = gca;
    set(h, "XTickLabel", [])
    set(h, "YTickLabel", [])
    xlabel("Posterior                                                                       Anterior")   % corresponding roughly to OpenSim X axis (horizontal pointing forward)
    ylabel("Inferior                                                                       Superior")      % corresponding to OpenSim Y axis (vertical pointing upwards)
    colorTitleHandle = get(hcb,'Title');
    titleString = 'time [s]';
    set(colorTitleHandle ,'String',titleString);

    %Plot Vector field of selected points
    quiver(-norm_fv_rotated(1:10:end,3),-norm_fv_rotated(1:10:end,1),-norm_fv_rotated(1:10:end,3), ...
    -norm_fv_rotated(1:10:end,1),"filled","AutoScale","on","Alignment","head", ...
   "ColorMode","auto","Marker","none","AutoScaleFactor",1,"LineWidth",1,"LineStyleMode","auto", ...
   "LineStyleMode","auto","MarkerEdgeColor","y","MarkerSize",10,"MaxHeadSize",10,"Color",	[0 0 0])

    hold off
    name_fig6svg = append(experiment_name, '_RCoPGH.png');
    name_fig6 = append(experiment_name, '_RCoPGH.svg');
    saveas(f6, name_fig6)
    saveas(f6, name_fig6svg)

    % plot the position of the GH force on the glenoid
    radius = sind(LmaxAngle);
    p=nsidedpoly(1000, 'Center', [0,0], 'Radius', radius);
    c = linspace(0,timesExp(end),length(Lnorm_fv_rotated));
    f7 = figure;
    hold on
    plot(p, 'FaceColor', 'r')
    for time_instant=1:numTimePoints
        scatter(-Lnorm_fv_rotated(time_instant,3), -Lnorm_fv_rotated(time_instant,1), [], c(time_instant), 'filled')
    end
    hcb = colorbar;
    h = gca;
    set(h, "XTickLabel", [])
    set(h, "YTickLabel", [])
    xlabel("Anterior                                                                       Posterior")   % corresponding roughly to OpenSim X axis (horizontal pointing forward)
    ylabel("Superior                                                                       Inferior")      % corresponding to OpenSim Y axis (vertical pointing upwards)
    colorTitleHandle = get(hcb,'Title');
    titleString = 'time [s]';
    set(colorTitleHandle ,'String',titleString);

    %Plot Vector field of selected points
    quiver(-Lnorm_fv_rotated(1:10:end,3),-Lnorm_fv_rotated(1:10:end,1),-Lnorm_fv_rotated(1:10:end,3), ...
    -Lnorm_fv_rotated(1:10:end,1),"filled","AutoScale","on","Alignment","head", ...
   "ColorMode","auto","Marker","none","AutoScaleFactor",1,"LineWidth",1,"LineStyleMode","auto", ...
   "LineStyleMode","auto","MarkerEdgeColor","y","MarkerSize",10,"MaxHeadSize",10,"Color",	[0 0 0])

    hold off
    name_fig7 = append(experiment_name, '_LCoPGH.png');
    name_fig7svg = append(experiment_name, '_LCoPGH.svg');
    saveas(f7, name_fig7)
    saveas(f7, name_fig7svg)

    %Plot the ratio of ligament length to ligament resting length
    f8=figure;
    plot(timesExp,lig_ratio,LineWidth=2)
    legend(ligaments,'Orientation','vertical','Location','best','Interpreter','none')
    title("Length Ratio")
    xlabel("Time [s]")
    ylabel("Ligament Length / Resting Length")
    saveas(f8, "Ligament_Length_Ratio.png")
    saveas(f8, "Ligament_Length_Ratio.svg")

    % f9=figure;
    % plot(timesExp,tb_Optim,LineWidth=2)
    % hold on
    % plot(timesExp,tA_Optim,LineWidth=2)
    % legend("Time to Fmincon", "Time Fmincon")
    % saveas(f9, "Time_Fmincon.png")
    % saveas(f9, "Time_Fmincon.svg")
end

%% SAVING THE RESULTS TO FILE
name_file = append('muscle_activations_', subject_considered, '_', experiment_name);

muscleNames = ArrayStr();
muscles.getNames(muscleNames);

muscle_order = "";
for i = 1:numMuscles
    muscle_order= [muscle_order, string(muscleNames.get(i-1).toCharArray')];
end

for i=1:length(coordNames)
    muscle_order= [muscle_order, string(coordNames{i})];
end

muscle_order = muscle_order(2:end);

% rescale the frequency of the solution knowing the freq of the data
frequency_solution = frequency_trc_data/time_interval;

% save(name_file,'-append')

save(name_file, 'xsol', 'muscle_order', 'frequency_solution', 'optimizationStatus', 'unfeasibility_flags', 'tOptim', ...
    'norm_fv_rotated', 'Lnorm_fv_rotated', 'A_eq_force', 'LA_eq_force', 'A_eq_acc', ...
    'coordinates','coordNames','speeds','accelerations','actsNames', ...
    'simulatedAccelerations','maxAngle', 'LmaxAngle', 'timesExp','numTimePoints', ...
    'lig_force','ligaments','lig_ratio','tb_Optim','tA_Optim','wc','force_vecc','Lforce_vecc');

file_results = fullfile(saving_path, name_file, '.mat');

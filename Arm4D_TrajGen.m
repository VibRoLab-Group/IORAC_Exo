%% Trajectory Generations for the Stationary Upper Limb Exoskeleton (Arm4D) simulations
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script randomly generates the tracking references and disturbances for all   
% stationary upper limb exoskeleton simulations

PathSetup; % Include folders to path directory and reset workspace

%% Trajectory Generations

% Remember to generate new trajectories when the simulation and control
% sample rates are changed!

Arm4D_TestNum=4; % This is only use to load Arm4D_Parameters
Arm4D_Parameters; % Load simulation parameter

% Use "help genFreqTraj" to view instructions.

%% For Link4 uncertainty only, no disturbance simulation (Fig.2)

tSpan=0:hh:300; % Reference time span

% Generate reference (position, velocity and acceleration) (rad);
[refData,refdtData,refddtData]=genFreqTraj(4,tSpan,0.4,linspace(0.2,0.4,10));

% Adjust the motion range of the third joint to be smaller and add
% constant joint initial position
refData=refData.*[1;1;0.5;1]+initialPositionVal*setInitialPosition;
refdtData=refdtData.*[1;1;0.5;1]; 
refddtData=refddtData.*[1;1;0.5;1];

% Save reference to file
save(strcat('./Data/simLink4Ref'),'refData','refdtData','refddtData');

%% For Link 3 and Link 4 uncertainties, no disturbance simulation (Fig.3)

tSpan=0:hh:600; % Reference time span

% Generate reference (position, velocity and acceleration) (rad);
[refData,refdtData,refddtData]=genFreqTraj(4,tSpan,0.4,linspace(0.2,0.4,10));

% Adjust the motion range of the third joint to be smaller and add
% constant joint initial position
refData=refData.*[1;1;0.5;1]+initialPositionVal*setInitialPosition;
refdtData=refdtData.*[1;1;0.5;1];
refddtData=refddtData.*[1;1;0.5;1];

% Save reference to file
save(strcat('./Data/simLink34Ref'),'refData','refdtData','refddtData');

%% For Link 3 and Link 4 uncertainties, disturbance simulation (Fig.4)

tSpan=0:hh:300; % Reference time span
% Generate reference (position, velocity and acceleration) (rad);
[refData,refdtData,refddtData]=genFreqTraj(4,tSpan,0.4,linspace(0.2,0.3,10));

% Adjust the motion range of the third joint to be smaller and add
% constant joint initial position
refData=refData.*[1;1;0.5;1]+initialPositionVal*setInitialPosition;
refdtData=refdtData.*[1;1;0.5;1];
refddtData=refddtData.*[1;1;0.5;1];

% Generation of disturbance, noise, and tremor
% Important: extremely large disturbance, noise, and tremor requires
% increasing the simulation and control sample rates in parameters. 
% Otherwise, the simulation will crash due to numerical error.
[disturbData]=genFreqTraj(4,tSpan,1,linspace(4,4*sqrt(5),20)); %Disturbance Torque (N-m)

% Save reference to file
save(strcat('./Data/simLink34RefDisturbance'),'refData','refdtData','refddtData','disturbData');
%% Trajectory Generations for the Wearable Wrist Exoskeleton (TAWE) simulations
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script generates the reference, disturbance, noise, and tremor
% trajectories used for all wrist exoskeleton simulations

PathSetup; % Include folders to path directory and reset workspace

%% Trajectory Generations

% Remember to generate new trajectories when the simulation and control
% sample rates are changed!

TAWE_Parameters; % Load tSpan from simulation parameter

% Use "help genFreqTraj" to view instructions.

% Generate reference (position, velocity and acceleration) (rad);
[refData,refdtData,refddtData]=genFreqTraj(2,tSpan,0.4,linspace(0.2,0.4,10));
% Scale the wrist radial ulnar deviation reference to a reasonable wrist workspace range 
refData=refData.*[0.5;1]; 
refdtData=refdtData.*[0.5;1];
refddtData=refddtData.*[0.5;1];

% Generation of disturbance, noise, and tremor
% Important: extremely large disturbance, noise, and tremor requires
% increasing the simulation and control sample rates in parameters. 
% Otherwise, the simulation will crash due to numerical error.
tremorData=genFreqTraj(2,tSpan,2e-1,linspace(3,6,7)); % Tremor Excitation Torque (N-m)
disturbData=genFreqTraj(2,tSpan,1e-2,linspace(10,10*sqrt(5),20)); % Disturbance Torque (N-m)
noiseData1=genFreqTraj(3,tSpan,5e-4,linspace(20,20*sqrt(5),20)); % Displacement Measurement Noise (m)
noiseData2=genFreqTraj(3,tSpan,2e-3,linspace(20,20*sqrt(5),20)); % Euler Angle Measurement Noise (rad)
noiseData=[noiseData1;noiseData2];

% Save reference to file
save(strcat('./Data/simTAWERef'),'refData','refdtData','refddtData','tremorData','disturbData','noiseData');

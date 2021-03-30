%% Parameters for the Wearable Wrist Exoskeleton (TAWE) simulations
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script contains the simulation, model, and control parameters
% for all wearable wrist exoskeleton (TAWE) simulations. 


%% Simulation Parameters
simFreq=500; % Simulation sample rate
ctrlFreq=simFreq/2; % Control and WKI algorithm sample rate
%Rememeber to adjust sample rates if you would like to use larger control
%gains or tremor/disturbance/noise amplitudes/frequencies

hh=1/simFreq; % Simulation step size
tt=0; % Time variables
flowNum=1; % This is used to select different sets of model constraints, set to 1 for all TAWE simulations
tSpan=0:hh:150; % Default time span

%% Model Parameters

% System parameters
rTetra=0.2/3*sqrt(6); % Dimension of the load tetrahedron (for center of mass estimation) on the hand

% The parameters of P_Sys are: 
% (1) Gravity acceleration
% (2) Gravity acceleration of the uncertain body (set to 1 for simulation, 
%     set to 0 when load tetrahedron (for center of mass estimation) is used)
% (3) Soft constraint Stiffness (used to prevent accumulated drifting error
%     during numerical integration of constrained system. This term will
%     not be active if the constraint is not significantly violated)
% (4) Hand damping (Assumed to be unknown for controller design)
% (5) Exoskeleton damping (Assumed to be unknown for controller design)
% (6) Load tetrahedron (for center of mass estimation) dimension parameter.

P_Sys_Val=  [9.8067;  9.8067;  100;  0.1;  0.05;  rTetra]; % This is used in the simulation 
P_Sys_Ctrl= [9.8067;  0;       100;  0.0;  0.0;  rTetra]; % This is used for the controller design


% Kinematic parameters obtained from CAD designs, which matches with the 3D
% models in visualization. Modifications are not recommended.

% Constant angles in the forearm/wrist model
P_Init1_Val=1*pi*[0;0;0; -10;15;5; 0;0;-10]/180; 
% For controller  with WKI, the initial frame of the wrist joint is unknown
P_Init1_Ctrl=[P_Init1_Val(1:3);zeros(3,1);P_Init1_Val(7:9)]; 

% Constant angles/offsets in TAWE (known)
P_Init2_Val=[0;0;0;-0.3699;0.7732;-1.8355;0.7874;0.1527;0.3268];
P_Init2_Ctrl=P_Init2_Val;

% Translational displacements in the wrist model obtained from CAD
P_Arm_Val=[-6.28; 70; -41; ...
           11.76; -60.87; -27.98;...
           -5; -30; 1]*1e-3; 
P_Arm_Val(4:6)=-ypr2Mat(P_Init1_Val(7:9))*P_Arm_Val(4:6); % A parameter adjustment to match with 3D model
% For controller with WKI, the first two 3D displacements are unknown
P_Arm_Ctrl=[zeros(6,1);...
           -5; -30; 1]*1e-3; 
% P_Arm_Ctrl(4:6)=-ypr2Mat(P_Init1_Ctrl(7:9))*P_Arm_Ctrl(4:6);

% TAWE mechanism link parameters obtained from CAD (known)
P_Exo_Val=[-68.6; 12; 21.5; ...
           120;120;30;10;-25]*1e-3;
P_Exo_Ctrl=P_Exo_Val;

% Hand Center of Mass displcaement parameters obtained from CAD
P_COM1_Val=[
            -15.38; 16.49; -28.05
            ]*1e-3;
% Assumed to be unknown for controller
P_COM1_Ctrl=zeros(3,1);
        
% TAWE mechanism Center of Mass displcaement parameters obtained from CAD (known)
P_COM2_Val=[
             -3.08; -10.54; -1.87;
            20.42; 48.05; 0;
            34.56; 60; 0;
            9.06; 8.28; 0;
            1.72; -4.49; -10.65;
            -1.32; -18.35; -2.32;
            ]*1e-3;
P_COM2_Ctrl=P_COM2_Val;
        
% Hand mass and moment of inertia obtained from CAD
% columns: (1) Mass; (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
P_Inert1=[
          0.656;    0.656*[+2.259e-03;-5.823e-04;+1.027e-04;+1.430e-03;+2.206e-04;+3.354e-03];
          ];
% Assumed to be unknown for controller
P_Inert1_Ctrl=zeros(7,1);
      
% TAWE mechanism link masses and moments of inertia obtained from CAD (known)
% columns: (1) Mass; (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
P_Inert2=[
          0.092;    0.092*[+2.267e-04;-1.115e-05;-3.116e-06;+1.831e-04;-9.881e-06;+2.478e-04];  ...
          0.028;    0.028*[+1.490e-03;-1.936e-04;-2.869e-04;+1.114e-03;-6.778e-04;+5.672e-04];  ...
          0.032;    0.032*[+2.169e-03;-1.516e-04;+3.318e-04;+1.860e-03;+8.138e-04;+4.813e-04];  ...
          0.012; 	0.012*[+5.682e-05;+2.011e-05;+9.749e-09;+1.630e-04;-2.356e-09;+1.776e-04];  ...
          0.018;    0.018*[+1.602e-04;+4.302e-07;-4.109e-05;+1.493e-04;-2.381e-05;+1.656e-04];  ...
          0.030;    0.030*[+2.260e-04;+1.539e-05;+1.650e-12;+3.106e-04;-6.426e-12;+5.167e-04];  ...
         ];
P_Inert2_Ctrl=P_Inert2;

% Uncertain load mass and inertia (set to zeros for simulation)
P_Uncertain_Val=zeros(7,1);
P_Uncertain_Ctrl=zeros(7,1);
      
% All parameters used for simulation, controller without WKI, and
% controller with WKI
ppValSim=[P_Sys_Val;P_Arm_Val;P_Exo_Val;P_Init1_Val;P_Init2_Val;P_COM1_Val;P_COM2_Val;P_Inert1;P_Inert2;P_Uncertain_Val];
ppValCtrlNoWKI=[P_Sys_Ctrl;P_Arm_Val;P_Exo_Val;P_Init1_Val;P_Init2_Val;P_COM1_Val;P_COM2_Val;P_Inert1_Ctrl;P_Inert2_Ctrl;P_Uncertain_Ctrl]; % This set assumes we know the kinematics
ppValCtrlWKI=[P_Sys_Ctrl;P_Arm_Ctrl;P_Exo_Ctrl;P_Init1_Ctrl;P_Init2_Ctrl;P_COM1_Ctrl;P_COM2_Ctrl;P_Inert1_Ctrl;P_Inert2_Ctrl;P_Uncertain_Ctrl];

ppNum=numel(ppValCtrlWKI); % Total model parameter count
ppCtrlCertainRange=[1:6,13:26,30:ppNum-7]; % Known parameter No. for controller design
ppKineUncertainRange=[7:12,27:29]; % Uncertain wrist kinematic Paramter No, to be estimated by WKI
ppDynaUncertainRange=[ppNum-6:ppNum]; % Uncertain inertia parameter No. to be estimated by adaptive control
% ppBaseOrientationRange=[24:26]; 

%% Control Parameters

% Adjust c1 for the disturbance attentuation gain, adjust K1 for the ratio
% between P and D gains, and adjust K3 for the value of K2.
c1Coeff=2; % c1 in IO-RAC
KK1Gain=4*eye(2); % K1 in IO-RAC
KK3Gain=2*eye(2); % K3 from IO-RAC Lyapunov function 
KK4Gain=0.5*KK3Gain/KK1Gain; % K4 from IO-RAC Lyapunov function
KK2Gain=2*KK4Gain/c1Coeff; % K2 in IO-RAC, obtained through stability proof

% BMFLC reference model parameters
bmflcRes=16; % Resolution, i.e., frequency components
bmflcRange=(linspace(2,7,bmflcRes)*2*pi).'; % Frequencies from Bandwidth (2-7 Hz by default)

% Adaptive parameter update gains
GammaInv1=diag([1 ones(1,6)*1e-1 ones(1,3)]); % Update Gain for uncertain inertia and load
GammaInv2=eye(9)*1e1; % Update Gain for a damping compensator
GammaInv3=eye(bmflcRes*4)*1e1; % Update Gain for BMFLC amplitudes


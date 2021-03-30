%% Parameters for the Stationary Upper Limb Exoskeleton (Arm4D) simulations
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script contains the simulation, model, and control parameters
% for all Stationary Upper Limb (Arm4D) simulations. 


%% Simulation Parameters

simFreq=500; % Simulation sample rate
ctrlFreq=simFreq/2; % Control and WKI algorithm sample rate
%Rememeber to adjust sample rates if you would like to use larger control
%gains or tremor/disturbance/noise amplitudes/frequencies

hh=1/simFreq; % Simulation step size
tt=0; % Time variables
flowNum=1; % This is used to select different sets of model constraints, set to 1 for all TAWE simulations
tSpan=0:hh:300; % Default time span

setInitialPosition=1; % set to 1 for initial Joint  poisiton
initialPositionVal=[pi/6;-pi/3;-pi/4;pi/3]*setInitialPosition; % initial Joint potion 

%% Model Parameter

% The following parameters are obtained from CAD
% Base to joint 1 displacement
link0=[50; (150+50*sqrt(3)); 0]*1e-3/2;
  
link1Inertia=[78.573*1e6;7939685;-11098819;-223;26236314;-58;33678890;]*1e-6/4; % Link 1 inertia (1) Mass; (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
linkCOM1=[-510.725;219.984;-0.011;]*1e-3/2; % Link 1 center of mass (COM) displacement 
link1=[-550+50;650;0;]*1e-3/2; % Joint 1 to Joint 2 displacement (dimension of Link 1)
  
link2Inertia=[15.512*1e6;731721;-42753;35961;199715;70172;823340;]*1e-6/4; % Link 2 inertia (1) Mass; (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
linkCOM2=[-29.124;170.018;14.494;]*1e-3/2; % Link 2 center of mass (COM) displacement 
link2=[150;320;0;]*1e-3/2; % Joint 2 to Joint 3 displacement (dimension of Link 2)

link3Inertia=[25279*1e6;1567217739;-456248765;5711433;447860306;-11489155;1834530423;]*1e-9/4; % Link 3 inertia (1) Mass; (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
linkCOM3=[-69.24;192.48;3.03;]*1e-3/2; % Link 3 center of mass (COM) displacement 
link3=[-150;420;0;]*1e-3/2; % Joint 3 to Joint 4 displacement (dimension of Link 3)
  
link4Inertia=[20702*1e6;3610423705;1310018301;-14622477;723995582;-22952557;4246343915;]*1e-9/4; % Link 4 inertia (1) Mass; (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
linkCOM4=[164.86;324.08;-5.32;]*1e-3/2; % Link 4 center of mass (COM) displacement 
link4=[0;400;0;]*1e-3/2; % Joint 4 to Exoskeleton end-point displacement (dimension of Link 4)
  
  
% Estimated Link 3 Inertia and COM (without the human body)
link3InertiaEst=[10.095*1e6;859486;-436843;1618;352967;1773;1143175;]*1e-6/4;
linkCOM3Est=[-156.631;232.450;15.714;]*1e-3/2;
link3InertiaEst([3 4 6])=link3Inertia([3 4 6]); % Simplification by only consider the diagonal components

% Estimated Link 4 Inertia and COM (without the human body)
link4InertiaEst=[2.684*1e6;121101;6238;0;2134;0;120648;]*1e-6/4;
linkCOM4Est=[15.352;154.956;0.000;]*1e-3/2;
link4InertiaEst([3 4 6])=link4Inertia([3 4 6]); % Simplification by only consider the diagonal components
  

% For the first simulation, only assume link 4 has unknown body
if Arm4D_TestNum==4
    link3InertiaEst=link3Inertia; %Activate for 1-body test
    linkCOM3Est=linkCOM3; %Activate for 1-body test
end
     
ppInertiaValSim=[link1Inertia;link2Inertia;link3Inertia;link4Inertia]; %Inertia parameters used for simulation
ppInertiaValCtrl=[link1Inertia;link2Inertia;link3InertiaEst;link4InertiaEst]; %Inertia parameters used for control
ppDimValSim=[link0;link1;link2;link3;link4;linkCOM1;linkCOM2;linkCOM3;linkCOM4]; %Kinematic parameters used for simulation
ppDimValCtrl=[link0;link1;link2;link3;link4;linkCOM1;linkCOM2;linkCOM3Est;linkCOM4Est]; %Kinematic parameters used for control
  
ang0=[0;0;0]; % Base orientation (Euler Angle), assumed to be zeros
ppAngleVal=[ang0];

% Unknown Inertia parameters 
% (1-4 for Link 3 uncertain mass and moments (xx,yy,zz), 
% 5-8 for Link 4 uncertain mass and moments (xx,yy,zz))
% Designed as defined in the paper.
ppInertiaUnknownVal=zeros(8,1); 
% Unknown mass parameters for load tetrahedron (for estimation of COM)
% (1-3 for Link 3, 4-6 for Link 4)
ppLoadUnknwonVal=zeros(6,1);
% The uncertain parameters of link 3 will not be effective for the first
% simulation where only link 4 has uncertain body

gAcc=9.81*1;
cJoint=1e-2;
rTetra=0.25/3*sqrt(6);
ppMiscVal=[gAcc;cJoint;rTetra]; % Gravitational acceleration, joint damping, and load tetrahedron dimension

% Model parameters for simulation
ppSimVal=[ppInertiaValSim;ppDimValSim;ppAngleVal;ppInertiaUnknownVal;ppMiscVal];
ppCtrlVal1=[ppInertiaValCtrl;ppDimValCtrl;ppAngleVal]; % Model paramters for controller desing (Part 1)
ppCtrlVal2=[ppInertiaUnknownVal]; % Model paramters for controller desing (Part 2)
ppCtrlVal3=[ppMiscVal]; % Model paramters for controller desing (Part 3)

%% Control Parameters

% Adjust c1 for the disturbance attentuation gain, adjust K1 for the ratio
% between P and D gains, and adjust K3 for the value of K2.
% Parameters match with Table 1 in the paper
c1Coeff=2; % c1 in IO-RAC
KK1Gain=4*eye(4); % K1 in IO-RAC
KK3Gain=2*eye(4); % K3 from IO-RAC Lyapunov function 
KK4Gain=0.5*KK3Gain/KK1Gain; % K4 from IO-RAC Lyapunov function
KK2Gain=2*KK4Gain/c1Coeff; % K2 in IO-RAC, obtained through stability proof

% Adaptivve parameter update gain
GammaInv=2*diag([1 ones(1,3)*1e-1 1 ones(1,3)*1e-1 ones(1,6)*1e0]);
if Arm4D_TestNum == 4
    GammaInv=2*diag([ 1 ones(1,3)*1e-1 ones(1,3)*1e0]);%Activate for link 4 uncertainy only test
end

% Values for checking inertia parameter estimation errors
CheckVal1=[link3Inertia([1 2 5 7])-link3InertiaEst([1 2 5 7]);link4Inertia([1 2 5 7])-link4InertiaEst([1 2 5 7])];

% Calculates the tetrahedron vertice translational positions
tretraVertices3 = [[-1;-1/sqrt(6);-1/sqrt(3)] [1;-1/sqrt(6);-1/sqrt(3)] [0;-1/sqrt(6);2/sqrt(3)] [0;3/sqrt(6);0]]*rTetra + linkCOM3Est;
tretraVertices4 = [[-1;-1/sqrt(6);-1/sqrt(3)] [1;-1/sqrt(6);-1/sqrt(3)] [0;-1/sqrt(6);2/sqrt(3)] [0;3/sqrt(6);0]]*rTetra + linkCOM4Est;
      
% Calculate the COM of the uncertain body, these are used to check if the
% COM estimation performance from Link 3 and Link 4 load tetrahedrons 
link3UncertainCOM = (link3Inertia(1)*linkCOM3 - link3InertiaEst(1)*linkCOM3Est)/(link3Inertia(1)- link3InertiaEst(1));
link4UncertainCOM = (link4Inertia(1)*linkCOM4 - link4InertiaEst(1)*linkCOM4Est)/(link4Inertia(1)- link4InertiaEst(1));

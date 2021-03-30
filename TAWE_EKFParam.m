%% EKF parameters for the Wearable Wrist Exoskeleton (TAWE) simulations
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script contains the default EKF parameters for the wrist kinematic
% identifications in the control simulation.


%% Default EKF parameters

% xx=[quat0;cxixi;dd0;dd1;kappa;kappaInt] where:
% (1) quat0 - initial alignment orientation of the base (socket) of wrist joint 
% (2) cxixi - a parameter in the wrist rotation constraint that determines represented the wrist model
% (3) dd0 - translational displacement from the forearm attachement (IMU1) to the base of wrist joint
% (4) dd1 - translational displacement from wrist joint to the hand attachement (IMU2)
% (5) kappa - the wrist joint rotation angles (respectively: Flexion-extension, 
%             Pronation-supination (constrained), and radial ulnar deviation)
% (6) kappaInt - integral of kappa, this term is used to increase the stability
%                of WKI, to make sure that estimated quat0 results in the mean 
%                values of kappa is at zeros, so as to prevent drifting of 
%                estimation

% pp=[stableMag;t;deltaT] where 
% (1) stableMag - set this term to zero to deactivate the effect of kappaInt
% (2) t - time variable
% (3) deltaT - EKF step size

% uu=[quat1;ddRef] where:
% (1) quat1 - rotation between IMU1 and IMU2
% (2) ddRef - translational displacement between IMU1 and IMU2

taweWKI.state.xx=double(0*taweWKI.info.xxSym); % State of the model (posteriori of the last step)
taweWKI.state.pp=double(0*taweWKI.info.ppSym); % Parameters of the model (will not be updated)
taweWKI.state.uu=double(0*taweWKI.info.uuSym); % Input of the model (e.g., sensor measurement, control input)

% Turn on the effects of kappaInt
taweWKI.state.pp(1)=1; 
taweWKI.state.pp(2)=1;
taweWKI.state.pp(3)=1;
taweWKI.state.pp(end)=1/ctrlFreq; % EKF time step is the same as control time step

taweWKI.state.xx(1)=1; % Real part of quat0 is set to 1 for zero-rotation
taweWKI.state.uu(1)=1; % Real part of quat1 quaternion is set to 1 for zero-rotation

% xxPlus=[quat0;cxixi;dd0;dd1;kappaExpr;kappaInt+kappaExpr*deltaT]
% This is the priori state prediction. Since estimated parameters are assumed
% constant or slowly time-varying, their priori prediction equals to their
% current values. 
% "kappaExpr" is an expression that calculates kappa based
% on quat0 and quat1: 
% kappaExpr=quat2YPR(quatMultiply(quatConj(quat0),quat1));
% Finally, kappaInt is obtained by integrating kappa along time.

% yy=[rotRegFull;transRegBase-ddRef;quatCons;quatStable];
% This is the priori observation (true values are zeros), where 
% (1) rotRegFull is the quaternion based wrist rotational constraint:
%     rotRegBase = 0 = [1 0 0 0]*quatMultiply([0;[0;1;0]],ypr2Quat(kappa))+cxixi*sin(kappa(1)/2)*sin(kappa(3)/2);
% (2) transRegBase - ddRef identifies the wrist translational displacement
%     transRegBase - ddRef = 0 = dd0+quat0Mat*kappaMat*dd1 - ddRef;
%     where quat0Mat and kappaMat are the rotation matrix of quat0 and kappa
% (3) quatStable = 0 = stableMag.*kappaInt = 0 is used to ensure that estimated 
%     quat0 results in the mean values of kappa is at zeros, so as to prevent 
%     drifting of estimation
taweWKI.state.yyHat=double(0*taweWKI.info.yySym); % observation based on priori State

taweWKI.state.PPxx=(1e-6*eye(numel(taweWKI.state.xx))); % Covariance matrix of the state (posteriori of the last step)
taweWKI.state.QQ=(1e-6*eye(numel(taweWKI.state.uu))); % Covariance of process noise
taweWKI.state.RR=(1e-6*eye(numel(taweWKI.state.uu))); % Covariance of observation noise

taweWKI.state.QQBase=(1e-10*eye(numel(taweWKI.state.xx))); % Auxiliary linear covariance of process noise to avoid singularity
taweWKI.state.QQBase(1:5,1:5)=taweWKI.state.QQBase(1:5,1:5)*1e2;
taweWKI.state.RRBase=(1e-6*eye(numel(taweWKI.state.yyHat))); % Auxiliary linear covariance of observation noise to avoid singularity
taweWKI.state.RRBase(6:end,6:end)=1e6*taweWKI.state.RRBase(6:end,6:end);
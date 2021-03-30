%% The Wearable Wrist Exoskeleton (TAWE) Control Simulation - without WKI
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script contains the simulation code for the wearable wrist
% exoskeleton with wrist kinematic identification and disturbance. The
% simulations also feature tremor and its suppression via the proposed controller.

PathSetup; % Include folders to path directory and reset workspace

%% 3D Visualization (Uncomment to load 3D models)
% Below is the initialization codes for the 3D visualization of the
% wrist exoskeleton, which includes 3D models, load tetrahedron (for center 
% of mass estimation), and main coordinate frames. Uncomment the code in this 
% section to initiate visualization. (Requires OpenGL for smooth animation)

%{
Sim=sAxes('TAWE Simulation',3,'Frames_TAWE.mat',@numTF_TAWE_mex);
Sim.setAxesProp('ArmIMU1Frame',1*[-0.2 -0.1 -0.2;0.2 0.5 0.2],[210 15]).setPresetTF('WORLD');

[ArmPose,ExoPose,HandCOMTraj,HandTetra]...
    =Sim.genPlot({'ArmPose';'ExoPose';'HandCOMTraj';'HandTetra'});
ArmPose.setLineSpec('-','none','r',2).setPlotPoint({'ArmIMU1Frame';'ArmDevFlexFrame';'ArmIMU2Frame'});
ExoPose.setLineSpec('-','none','b',2).setPlotPoint({'ExoIMU1Frame';'ExoMedFrame';'ExoPanFrame';'ExoLink12Frame';...
                                                    'ExoLink23Frame';'ExoLink34Frame';'ExoLink45Frame';'ExoAttachFrame';'ExoIMU2Frame'});
HandTetra.setLineSpec(':','+','g',3).setPlotPoint({'Handmpt1Frame';'Handmpt2Frame';'Handmpt3Frame';'Handmpt4Frame';'Handmpt2Frame';'Handmpt1Frame';'Handmpt3Frame';'Handmpt1Frame';'Handmpt4Frame';});

[ArmBaseCAD,ArmHandCAD...
     ,ExoBaseCAD,ExoPanCAD,ExoLink1CAD,ExoLink2CAD,...
     ExoLink3CAD,ExoLink4CAD,ExoLinkJointCAD...
    ]...
=Sim.genPatch({
               'ArmBaseCAD'  'ArmIMU1Frame';...
               'ArmHandCAD'  'ArmIMU2AdjustFrame';...
               'ExoBaseLink'  'ExoIMU1Frame';...
               'ExoPan'  'ExoPanFrame';...
               'ExoLink1'  'ExoLink12Frame';...
               'ExoLink2'  'ExoLink23Frame';...
               'ExoLink3'  'ExoLink34Frame';...
               'ExoLink4'  'ExoLink45Frame';...
               'ExoLinkJoint'  'ExoAttachFrame';...
               });
ArmBaseCAD.setFaceProp([0.8 0.5 0.5],0.5).setModel('taweArmBase.STL',0.1,[],0.001);
ArmHandCAD.setFaceProp([0.8 0.5 0.5],0.5).setModel('taweArmHand.STL',0.1,[],0.001);
ExoBaseCAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoBaseLink.STL',0.1,[],1);
ExoPanCAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoPan.STL',0.1,[],1);
ExoLink1CAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoLink1.STL',0.025,[],0.001);
ExoLink2CAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoLink2.STL',0.025,[],0.001);
ExoLink3CAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoLink3.STL',0.025,[],0.001);
ExoLink4CAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoLink4.STL',0.025,[],0.001);
ExoLinkJointCAD.setFaceProp([0.5 0.5 0.8],0.5).setModel('taweExoLinkJoint.STL',1,[],0.001);
Sim.Patch(1).setFaceProp([0.8 0.8 1],1).setModel([],1,[],0.1);
%}

%% Load Model Parameters

TAWE_Parameters; % This loads the TAWE simulation parameters

% This overwrites ppVal (specifically P_Init2 for the initial angles of the
% TAWE mechanism so that the closed kinematic chain is fulfilled)
% load taweCloseLoopParam.mat; 

% This loads the tracking reference (position, velocity, acceleration) and
% synthetic tremors, disturbances, and noises. The reference can be
% generated in TAWE_Simulation_Setup
load('simTAWERef','refData','refdtData','refddtData','tremorData','disturbData','noiseData');
noiseFlag=1; % Set to 0 to turn off noise
disturbFlag=1; % Set to 0 to turn off disturbance

%% Initialize WKI with EKF for sometime before the controller
% This is important since the all initial kinematic parameters are set to
% zeros, which will break the controller if directly applied to controller

tSpan1=0:hh:30; % Initialization Time Span

% State vector qq and qqdot, which are:
% (1) Wrist Radial Ulnar Deviation (RUD) (rad) Rotation on z axis in the wrist frame
% (2) Wrist Flexion Extension (FE) (rad) Rotation on x axis in the wrist frame
% (3) Wrist Pronation Supination (SP) (rad), which is constrained to FE and
%     RUD, Rotation on y axis in the wrist frame
% (4-6) TAWE Mechanism Joint Angles (rad)
simStateVal=zeros(18,1);

% The input vectors uu are defined as 
% (1-2) human inputs at the wrist FE and RUD directions
% (3-4) control input from TAWE at the first two links
% (5-8) virtual input at the load tetrahedron vertices (-z direction) 
%       (for center of mass estimation), used to acquire the Jacobian for 
%       uncertain load parameter more conveniently. These inputs are not 
%       used for simulation, but for controller design only.
simInputVal=zeros(8,1); % Used for simulation 

% Nonholonomic state used in constraint only. Note that are not internal 
% states of the system, but Euler Angles qq_nh that satisfies
% RotMat(qq1) = RotMat(qq_nh)*RotMat(qq_nh)
% Hence, the rotation by qq1 is split to two identical rotations by qq_nh
% in sequance. The wrist kineamtic model uses this to form a z-x-z-x 
% sequenced rotation joint, where all the z angles are the same, and all
% the x angles are the same. This assumes that the wrist rotation is
% devided evenly into the radiocarpal and midicarpal joints.
nhSignal=zeros(3,1);

% Constraint force for simulation process (not used in controller design)
consForce=zeros(7,1); 

load taweWKI.mat; % Load WKI model
TAWE_EKFParam; % Initialize WKI parameters

tic;
for ii=1:numel(tSpan1)
    % Here we use the references generated starting from t=2 min, so  that 
    % it does not overlap with the control simulation reference
    ll=ii+simFreq*120; 
    tt=tSpan1(ii);
    
    qq=simStateVal(1:end/2);
    qqdot=simStateVal(end/2+1:end);
    
    ref=refData(:,ii);
    refdt=refdtData(:,ii);
    refddt=refddtData(:,ii);
    
    % Here we assume that the user provides the control input, hence all
    % true model parameters are assumed to be available. The controller is
    % a simple model-based feedback controller. The control performance is
    % not important, as we only care about the generated position trajectories  
    % for wrist kinematic identification.
    
    % Get System Properties
    [~,~,~,~,framePos,frameQuat,MM,~,GFF,~,JacU,JacCons1,~,~,CenMat,InertLeftComponent,InertRightComponent1]...
        =System_TAWE_mex(tSpan1(ii),simStateVal,ppValSim,simInputVal,nhSignal);
    MM=sum(MM,3); % Sum up inertia matrix of all bodies
    JacU=sum(JacU,3).'; % Sum up constraint Jacobian Matrices of all inputs
    JacCons1=sum(JacCons1,3); % Sum up constraint Jacobian Matrices of all constraints
    GFF=sum(GFF,2); % Sum up generalized forces (excluding those from inputs and load tetrahedron (for center of mass estimation))
    CenMat=sum(CenMat,3); % Sum up all centripetal matrices'

    % Properties used to numerically calculate Coriolis matrix and 
    % constraint Jacobian derivative
    simStateVal2=simStateVal;
    % Use the current state velocity to estimate the position after a
    % short time step
    simStateVal2(1:end/2)=simStateVal2(1:end/2)+simStateVal(end/2+1:end)*(hh); 

    % Calculate the equivalent qq_nh used in the constraint (see
    % above), Note that are not internal states of the system, but Euler 
    % Angles qq_nh that satisfies
    % RotMat(qq1) = RotMat(qq_nh)*RotMat(qq_nh)
    % Hence, the rotation by qq1 is split to two identical rotations by qq_nh
    % in sequance.
    rotQuat=ypr2Quat(simStateVal2([2 3 1])); % Unit quaternion equivalent to wrist rotation
    halfCos=sqrt((rotQuat(1)+1)/2); 
    rotQuatSquareRoot=[halfCos;rotQuat(2:4)/2/halfCos]; % Half of the rotQuat rotation
    nhSignal2=quat2YPR(rotQuatSquareRoot)*2; % Calculate qq_nh estimate after a short time step
    % Obtain estimated properties after a short time step
    [~,~,~,~,~,~,~,~,~,~,~,JacCons2,~,~,~,~,InertRightComponent2]=...
        System_TAWE_mex(tSpan1(ii),simStateVal2,ppValSim,simInputVal,nhSignal2);
    
    % Formulate constraint properties (please also see supplemental document)
    Jlam1=JacCons1(:,1:2); % The Jacobian for preseved coordinates of the constrained model 
    Jlam2=JacCons1(:,3:end); % The Jacobian for internal coordinates of the constrained model 
    JlamMap=-Jlam2\Jlam1; % The mapper between preseved coordinates and internal coordinates
    JlamMapFull1=[eye(2);JlamMap]; % The mapper between preseved coordinates and all coordinates

    % This step is similar to above, except using the estimated
    % properties
    Jlam1=JacCons2(:,1:2);
    Jlam2=JacCons2(:,3:end);
    JlamMap=-Jlam2\Jlam1;
    JlamMapFull2=[eye(2);JlamMap];

    % Numerically calculate derivatives of matrices
    JlamMapFull=JlamMapFull1;
    JlamMapFullDot=(JlamMapFull2-JlamMapFull1)/(hh);
    InertRightComponent=InertRightComponent1;
    InertRightComponentDot=(InertRightComponent2-InertRightComponent1)/(hh);
    
    % Numerically calculate coriolis matrix
    CorMat=InertLeftComponent*InertRightComponentDot;
        
    % Calculated properties of constrained system
    MMRightCombine = JlamMapFull.'*MM;
    MMCombine = MMRightCombine*JlamMapFull; % Inertia matrix 
    CCCombine = JlamMapFull.'*(CorMat+CenMat)*JlamMapFull ...
                + MMRightCombine*JlamMapFullDot; % Coriolis and Centripetal Matrix
    GFFCombine = JlamMapFull.'*GFF; % Generalized force (which excludes estimated loads)
    JacU1Tpose = JlamMapFull.'*JacU(1:2,:).'; % Jacobian of user input (transposed)

    %Calculate the model based PD controller
    uPDterm=(ref-qq(1:2))+(refdt-qqdot(1:2));
    userInput=(JacU1Tpose)\(MMCombine*refddt+MMCombine*uPDterm...
                            + CCCombine*qqdot(1:2) - GFFCombine);
    simInputVal(1:2)=userInput;
    

    if rem(ii,simFreq/ctrlFreq)==0
        
        imu1Pos=framePos(:,2); % Get position measurement from IMU1 on the forearm
        imu2Pos=framePos(:,5); % Get position measurement from IMU2 on the hand dorsum
        imu1Quat=frameQuat(:,2); % Get rotation quaternion measurement from IMU1 on the forearm
        imu2Quat=frameQuat(:,5); % Get rotation quaternion measurement from IMU2 on the hand dorsum
        
        % Calculate displacement between frames with noise added
        posDiff = quat2Mat(imu1Quat).' * (imu2Pos-imu1Pos) + noiseFlag*noiseData(1:3,ll); 
        % Calculate rotation between frames with noise added
        rotDiff = quatMultiply(quatConj(imu1Quat),imu2Quat);
        rotDiff = ypr2Quat(quat2YPR(rotDiff) + noiseFlag*noiseData(4:6,ll));

        % Use sensor measurements as EKF input
        taweWKI.state.uu=[rotDiff;posDiff];
        % Iterate EKF
        taweWKI.state=ekfStep(taweWKI.func,taweWKI.state);
    end
    
    % Using Runge Kutta 4th method to integrate system ODE
    [simStateVal,nhSignal,consForce]=rk4(@Flow_TAWE_mex,hh,flowNum,tt,simStateVal,ppValSim,simInputVal,nhSignal,consForce);
end
toc;

taweWKI00=taweWKI; %Initial kinematic approximation after the WKI runs for 30 seconds



%% Simulation Loop (Uncomment to run)
% %{

tSpan2=0:hh:120;
simStateData=zeros(9,numel(tSpan2),3);

tic;


% Three different simulation loops:
% (jj=1) with WKI and disturbance
% (jj=2) with WKI, disturbance, and tremor
% (jj=3) with WKI, disturbance, tremor, and BMFLC model for tremorsuppression
for jj=1:3
    
    taweWKI=taweWKI00; % Load initialized WKI parameters
    
    % State vector qq and qqdot, which are:
    % (1) Wrist Radial Ulnar Deviation (RUD) (rad) Rotation on z axis in the wrist frame
    % (2) Wrist Flexion Extension (FE) (rad) Rotation on x axis in the wrist frame
    % (3) Wrist Pronation Supination (SP) (rad), which is constrained to FE and
    %     RUD, Rotation on y axis in the wrist frame
    % (4-6) TAWE Mechanism Joint Angles (rad)
    simStateVal=zeros(18,1); 
    ctrlStateVal=zeros(18,1); % Used for controller design since the wrist motion is estimated
    
    simParamVal=ppValSim; % Simulation uses true parameters
    ctrlParamVal=ppValCtrlWKI; % Controller uses estimated parameters
    
    simStateDataNow=zeros(9,numel(tSpan2)); % Used to collect data for the current jj
    refDataNow=refData;
    refdtDataNow=refdtData;
    refddtDataNow=refddtData;
    
    % Uncertain inertia and load parameter introduced by unknown body on the
    % hand, which are:
    % (1) Mass
    % (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
    % (8-10) Mass parameters on the load tetrahedron (for center of mass estimation)
    unCertainParam1=zeros(10,1);
    % Uncertain parameters for damping compensator at qqdot (9 elements)
    unCertainParam2=zeros(9,1);
    % Uncertain parameters for BMFLC compensator at qqdot (9 elements)
    unCertainParam3=zeros(bmflcRes*4,1); 
    
    % uncertain amplitude parameter Jacobian for BMFLC
    BMFLCJac=zeros(bmflcRes*4,2);
    
    tremorFlag=0; % Set to 0 to deactivate tremor
    BMFLCFlag=0; % Set to 0 to deactivate BMFLC
    if jj==2
        tremorFlag=1;
    end
    if jj==3
        tremorFlag=1;
        BMFLCFlag=1;
    end

    simInputVal=zeros(8,1); % Used for simulation 
    ctrllnputVal=zeros(8,1); % Used for controller calculation only
    userInput=zeros(2,1); % The user input vector
    taweInput=zeros(2,1); % The user input vector
    
    consForce=zeros(7,1); % Constraint force for simulation process (not used in controller design)

    % Nonholonomic state used in constraint only. Note that are not internal 
    % states of the system, but Euler Angles qq_nh that satisfies
    % RotMat(qq1) = RotMat(qq_nh)*RotMat(qq_nh)
    % Hence, the rotation by qq1 is split to two identical rotations by qq_nh
    % in sequance.
    nhSignal=zeros(3,1);
    
    for ii=1:numel(tSpan2)

        tt=tSpan2(ii); % Update time variable
        simStateDataNow(:,ii)=simStateVal(1:end/2);

        qq=simStateVal(1:end/2);
        qqdot=simStateVal(end/2+1:end);

        ref=refDataNow(:,ii);
        refdt=refdtDataNow(:,ii);
        refddt=refddtDataNow(:,ii);

        if rem(ii,simFreq/ctrlFreq)==0            
            

            % %%Kalman Filter Update 
            
            % Acquire kinematic properties (use simulation state as acquired by the senor)
            [~,frameVel,~,~,framePos,frameQuat]=System_TAWE_mex(tSpan2(ii),simStateVal,simParamVal,simInputVal,nhSignal);

            imu1Pos=framePos(:,2); % Get position measurement from IMU1 on the forearm
            imu2Pos=framePos(:,5); % Get position measurement from IMU2 on the hand dorsum
            imu1Quat=frameQuat(:,2); % Get rotation quaternion measurement from IMU1 on the forearm
            imu2Quat=frameQuat(:,5); % Get rotation quaternion measurement from IMU2 on the hand dorsum

            % Calculate displacement between frames with noise added
            posDiff = quat2Mat(imu1Quat).' * (imu2Pos-imu1Pos) + noiseFlag*noiseData(1:3,ll); 
            % Calculate rotation between frames with noise added
            rotDiff = quatMultiply(quatConj(imu1Quat),imu2Quat);
            rotDiff = ypr2Quat(quat2YPR(rotDiff) + noiseFlag*noiseData(4:6,ll));

            % Use sensor measurements as EKF input
            taweWKI.state.uu=[rotDiff;posDiff];
            % Iterate EKF
            taweWKI.state=ekfStep(taweWKI.func,taweWKI.state);
 
            % %%Controller Update 
            
            % Exoskeleton Joint angles are directly obtained from
            % simulation data
            ctrlStateVal(4:9)=qq(4:9);
            ctrlStateVal(13:18)=qqdot(4:9);
            
            % Obtain estimated wrist motion
            qq1Est=taweWKI.state.xx(12:14);
            ctrlStateVal(1:3)=[qq1Est(3);qq1Est(1);qq1Est(2)];
            
            % Convert measured angular velocity to wrist motion velocity (as Euler angle rate)
            imu2AngVel=quat2Mat(imu2Quat).'*frameVel(4:6,5); % Obtain local angular velocity 
            
            % A matrix to convert local angular velocity to yaw-pitch-roll
            % Euler angle rates
            angVel2EulRateMat=[
                            1,  sin(qq1Est(1))*tan(qq1Est(2)),  cos(qq1Est(1))*tan(qq1Est(2));
                            0,  cos(qq1Est(1)),                 -sin(qq1Est(1));
                            0,  sin(qq1Est(1))/cos(qq1Est(2)),  cos(qq1Est(1))/cos(qq1Est(2));
                          ];
                      
            % Obtain estimated wrist motion velocity
            qq1dotEst=angVel2EulRateMat*imu2AngVel;
            ctrlStateVal(10:12)=[qq1dotEst(3);qq1dotEst(1);qq1dotEst(2)];
            
            % Insert WKI estimated wrist kinematic parameters
            ctrlParamVal(ppKineUncertainRange)=[taweWKI.state.xx(6:11);quat2YPR(taweWKI.state.xx(1:4))];
            % Insert estimated inertia parameters
            ctrlParamVal(end-6:end)=unCertainParam1(1:7); 
            % Insert estimated mass parameters on load tetrahedron 
            % (for center of mass estimation), notice that the last one is 
            % calculated as the estimated total mass minus the other mass parameters
            ctrllnputVal(5:8)=[unCertainParam1(8:10); unCertainParam1(1)-sum(unCertainParam1(8:10))];

            
            % Calculate the equivalent qq_nh used in the constraint. Note 
            % that are not internal states of the system, but Euler 
            % Angles qq_nh that satisfies
            % RotMat(qq1) = RotMat(qq_nh)*RotMat(qq_nh)
            % Hence, the rotation by qq1 is split to two identical rotations by qq_nh
            % in sequance.
            rotQuat=ypr2Quat(ctrlStateVal([2 3 1]));
            halfCos=sqrt((rotQuat(1)+1)/2);
            rotQuatSquareRoot=[halfCos;rotQuat(2:4)/2/halfCos];
            nhSignal1=quat2YPR(rotQuatSquareRoot)*2; % Note that the estimation nhSignal1 is different from nhSignal obtained from simulation
            
            
            % Get System Properties
            [~,~,~,~,~,Quat,MM,GFI,GFF,GFU,JacU,JacCons1,~,~,CenMat,InertLeftComponent,InertRightComponent1]...
                =System_TAWE_mex(tSpan2(ii),ctrlStateVal,ctrlParamVal,ctrllnputVal,nhSignal1);
            MM=sum(MM,3); % Sum up inertia matrix of all bodies
            JacU=sum(JacU,3).'; % Sum up constraint Jacobian Matrices of all inputs
            GFF=sum(GFF,2); % Sum up generalized forces (excluding those from inputs and load tetrahedron (for center of mass estimation))
            GFU=sum(GFU(:,1:4),2); % Sum up load tetrahedron generalized forces (equivalent to estimated generalized uncertainload)
            CenMat=sum(CenMat,3); % Sum up all centripetal matrices'
            % Substitute the model generated wrist rotational constraint
            % Jacobian with the one estimated via WKI
            Jac7=[taweWKI.info.rotJacFunc(taweWKI.state.pp(5),ctrlStateVal([2 3 1])),zeros(1,6)]; 
            JacCons1(7,:)=Jac7; % Sum up constraint Jacobian Matrices of all constraints

            
            % Properties used to numerically calculate Coriolis matrix and 
            % constraint Jacobian derivative
            ctrlStateVal2=ctrlStateVal;
            % Use the current state velocity to estimate the position after a
            % short time step
            ctrlStateVal2(1:end/2)=ctrlStateVal2(1:end/2)+ctrlStateVal(end/2+1:end)*(hh);
            % Similar to above, calculate the equivalent qq_nh used in the constraint
            rotQuat=ypr2Quat(ctrlStateVal2([2 3 1]));
            halfCos=sqrt((rotQuat(1)+1)/2);
            rotQuatSquareRoot=[halfCos;rotQuat(2:4)/2/halfCos];
            % Note that the estimation nhSignal2 (after a short time step) is different from nhSignal obtained from simulation
            nhSignal2=quat2YPR(rotQuatSquareRoot)*2; 
            % Obtain estimated properties after a short time step
            [~,~,~,~,~,~,~,~,~,~,~,JacCons2,~,~,~,~,InertRightComponent2]=....
                System_TAWE_mex(tSpan2(ii),ctrlStateVal2,ctrlParamVal,ctrllnputVal,nhSignal2);
            % Substitute the model generated wrist rotational constraint
            % Jacobian with the one estimated via WKI
            Jac7=[taweWKI.info.rotJacFunc(taweWKI.state.pp(5),ctrlStateVal2([2 3 1])),zeros(1,6)]; 
            JacCons2(7,:)=Jac7;

            
            % Formulate constraint properties (please also see supplemental document)
            Jlam1=JacCons1(:,1:2); % The Jacobian for preseved coordinates of the constrained model 
            Jlam2=JacCons1(:,3:end); % The Jacobian for internal coordinates of the constrained model 
            JlamMap=-Jlam2\Jlam1; % The mapper between preseved coordinates and internal coordinates
            JlamMapFull1=[eye(2);JlamMap]; % The mapper between preseved coordinates and all coordinates

            
            % This step is similar to above, except using the estimated
            % properties
            Jlam1=JacCons2(:,1:2);
            Jlam2=JacCons2(:,3:end);
            JlamMap=-Jlam2\Jlam1;
            JlamMapFull2=[eye(2);JlamMap];

            
            % Numerically calculate derivatives of matrices (the Jacobian
            % matrices do not involve velocity (i.e. qqdot))
            JlamMapFull=JlamMapFull1;
            JlamMapFullDot=(JlamMapFull2-JlamMapFull1)/(hh);
            InertRightComponent=InertRightComponent1;
            InertRightComponentDot=(InertRightComponent2-InertRightComponent1)/(hh);

            
            % Numerically calculate coriolis matrix
            CorMat=InertLeftComponent*InertRightComponentDot;

            
            % Calculated properties of constrained system
            MMRightCombine = JlamMapFull.'*MM;
            MMCombine = MMRightCombine*JlamMapFull; % Inertia matrix 
            CCCombine = JlamMapFull.'*(CorMat+CenMat)*JlamMapFull ...
                        + MMRightCombine*JlamMapFullDot; % Coriolis and Centripetal Matrix
            GFFCombine = JlamMapFull.'*GFF; % Generalized force (which excludes estimated loads)
            JacU1Tpose = JlamMapFull.'*JacU(1:2,:).'; % Jacobian of user input (transposed)
            JacU2Tpose = JlamMapFull.'*JacU(3:4,:).'; % Jacobian of TAWE input (transposed)
            JacU3Tpose = JlamMapFull.'*JacU(5:8,:).'; % Jacobian of load tetrahedron forces (transposed)
            loadRotMat = quat2Mat(Quat(:,5)); % Frame rotation of the uncertain body (estimated)

            
            % Intermeditate terms calculated using the estimated state
            midTerm1=refddt-KK1Gain*(ctrlStateVal(10:11)-refdt); % derivative of zeta in the paper (Eq.(15))
            midTerm2=refdt-KK1Gain*(ctrlStateVal(1:2)-ref); % zeta in the paper (Eq.(15))
            midTerm3=ctrlStateVal(10:11) - refdt + KK1Gain*(ctrlStateVal(1:2)-ref); % xi in the paper (Eq.(15))
        

            
            % Control Input using estimated model properties
            JJww=0.5*JacU2Tpose.'; % estimated disturbance jacobian matrix
            RRInv=(JJww.'*JJww+KK2Gain); % R from the paper (Eq.(15))
            
            % The controller design from Eq.(14-16) using estimated model parameters
            taweInput=(JacU2Tpose)\( MMCombine*midTerm1... % compensate acceleration force
                           + CCCombine*midTerm2 ... % compensate coriolis force
                           - GFFCombine... % compensate generalized force
                           - JlamMapFull.'*GFU... % compensate uncertain load (from load tetrahedron (for center of mass estimation))
                           + JlamMapFull.'*(unCertainParam2.*qqdot)... % compensate uncertain damping
                           + JlamMapFull.'*JacU(1:2,:).'*BMFLCJac.'*unCertainParam3... % BMFLC for tremor suppression (not active with BMFLCFlag==0)
                           - c1Coeff * RRInv * midTerm3); % feedback term (Eq.(15))
        
            
            %Prepare for Uncertain Parameter Jacobian Calculation
            uncertainInertRightComponent=InertRightComponent(1:6,:)*JlamMapFull;
            uncertainInertRightComponentDot=InertRightComponentDot(1:6,:)*JlamMapFull+InertRightComponent(1:6,:)*JlamMapFullDot;
            
            % Jacobian for parameter from uncertain body inertia and load
            paramJac1 = uncertainBodyParamJacobianFunc_TAWE_mex...
                        (ctrlStateVal(10:11),midTerm1,midTerm2,loadRotMat,...
                        uncertainInertRightComponent,uncertainInertRightComponentDot,(JacU3Tpose.'));
            % Jacobian for parameter from uncertain damping
            paramJac2 = diag(qqdot)*JlamMapFull;
            % Numerically construct the Jacobian matrix for BMFLC
            if BMFLCFlag
                sineArray=sin(bmflcRange);
                cosineArray=cos(bmflcRange);
                BMFLCJac(1:bmflcRes,1)=sineArray;
                BMFLCJac(bmflcRes+1:2*bmflcRes,1)=cosineArray;
                BMFLCJac(2*bmflcRes+1:3*bmflcRes,2)=sineArray;
                BMFLCJac(3*bmflcRes+1:4*bmflcRes,2)=cosineArray;
                paramJac3=BMFLCJac;
            else
                BMFLCJac = BMFLCJac*0; % Turn off BMFLC if not in the active suppression simulation loop
            end
            
            % Update uncertain parameters, note that (simFreq/ctrlFreq) 
            % ensures the step size is the same for different control sample 
            % rates
            unCertainParam1 = unCertainParam1 - 1 * hh * GammaInv1 * paramJac1 * midTerm3 * (simFreq/ctrlFreq);
            unCertainParam2 = unCertainParam2 - 1 * hh * GammaInv2 * paramJac2 * midTerm3 * (simFreq/ctrlFreq);
            if BMFLCFlag
                unCertainParam3 = unCertainParam3 - 1 * hh * GammaInv3 * paramJac3 * midTerm3 * (simFreq/ctrlFreq);
            else
                unCertainParam3 = unCertainParam3*0; % Turn off BMFLC if not in the active suppression simulation loop
            end

        end
        
        userInput=tremorFlag*tremorData(:,ii); % Tremor excitation from user
        simInputVal(1:2)=userInput; % Insert into input for simulation
        simInputVal(3:4)=taweInput+disturbFlag*disturbData(:,ii); % Insert into input for simulation (overlaid with disturbance)

        % Using Runge Kutta 4th method to integrate system ODE
        [simStateVal,nhSignal,consForce]=rk4(@Flow_TAWE_mex,hh,flowNum,tSpan2(ii),...
                                                simStateVal,simParamVal,simInputVal,nhSignal,consForce);
                                            
        % The codes below refreshes the 3D visualization for run-time animation
        % Uncomment to run.
%         if rem(ii,15)==0 %adjust the frame rate through the number
%             Sim.drawNow(tSpan2(ii),simStateVal,simParamVal,simInputVal,nhSignal); % requires 3D environment
%         end
    end
    simStateData(:,:,jj)=simStateDataNow;
end

toc;
save(strcat('./Data/simTAWEDataWKI'),'simStateData');
%}
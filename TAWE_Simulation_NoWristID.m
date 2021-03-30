%% The Wearable Wrist Exoskeleton (TAWE) Control Simulation - without WKI
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script contains the simulation code for the wearable wrist
% exoskeleton without wrist kinematic identification, tremor, or disturbance.

PathSetup; % Include folders to path directory and reset workspace

%% 3D Visualization (Uncomment to load 3D models)
% Below is the initialization codes for the 3D visualization of the
% wrist exoskeleton, which includes 3D models, load tetrahedron (for center 
% of mass estimation), and main coordinate frames. Uncomment the code in 
% this section to initiate visualization. (Requires OpenGL for smooth animation)

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

%% Load Model Parameters (Do not comment)

TAWE_Parameters; % This loads the TAWE simulation parameters

% This loads the tracking reference (position, velocity, acceleration). The
% reference can be generated in TAWE_Simulation_Setup
load('simTAWERef','refData','refdtData','refddtData');


%% Simulation Initialization

tSpan1 = 0:hh:60; %Simulation Time Span

% State vector qq and qqdot, which are:
% (1) Wrist Radial Ulnar Deviation (RUD) (rad) Rotation on z axis in the wrist frame
% (2) Wrist Flexion Extension (FE) (rad) Rotation on x axis in the wrist frame
% (3) Wrist Pronation Supination (SP) (rad), which is constrained to FE and
%     RUD, Rotation on y axis in the wrist frame
% (4-6) TAWE Mechanism Joint Angles (rad)
simStateVal=zeros(18,1); 
simParamVal=ppValSim; % Simulation uses true parameters
ctrlParamVal=ppValCtrlNoWKI; % Controller uses estimated parameters

% Uncertain inertia and load parameter introduced by unknown body on the
% hand, which are:
% (1) Mass
% (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
% (8-10) Mass parameters on the load tetrahedron (for center of mass estimation)
unCertainParam1=zeros(10,1);
% Uncertain parameters for damping compensator at qqdot (9 elements)
unCertainParam2=zeros(9,1);


% The input vectors uu are defined as 
% (1-2) human inputs at the wrist FE and RUD directions
% (3-4) control input from TAWE at the first two links
% (5-8) virtual input at the load tetrahedron vertices (-z direction) 
%       (for center of mass estimation), used to acquire the Jacobian for 
%       uncertain load parameter more conveniently. These inputs are not 
%       used for simulation, but for controller design only.

simInputVal=zeros(8,1); % Used for simulation 
ctrllnputVal=zeros(8,1); % Used for controller calculation only
taweInput=zeros(2,1); % The TAWE input term

consForce=zeros(7,1); % Constraint force for simulation process (not used in controller design)

% Nonholonomic state used in constraint only. Note that are not internal 
% states of the system, but Euler Angles qq_nh that satisfies
% RotMat(qq1) = RotMat(qq_nh)*RotMat(qq_nh)
% Hence, the rotation by qq1 is split to two identical rotations by qq_nh
% in sequance. The wrist kinematic model uses this to form a z-x-z-x 
% sequenced rotation joint, where all the z angles (on RUD direction) are the same, and all
% the x angles (on FE direction) are the same. This assumes that the wrist rotation is
% devided evenly into the radiocarpal and midicarpal joints.
nhSignal=zeros(3,1);

simStateData=zeros(9,numel(tSpan1)); % Store state data from simulation



%% Simulation Loop
tic;
for ii=1:numel(tSpan1)
    
    tt=tSpan1(ii); % Update time variable
    simStateData(:,ii)=simStateVal(1:end/2); % Store qq
    
    qq=simStateVal(1:end/2);
    qqdot=simStateVal(end/2+1:end);
    
    ref=refData(:,ii);
    refdt=refdtData(:,ii);
    refddt=refddtData(:,ii);
    
    if rem(ii,simFreq/ctrlFreq)==0 % used for control sample rate
        % Insert estimated inertia parameters
        ctrlParamVal(end-6:end)=unCertainParam1(1:7); 
        % Insert estimated mass parameters on load tetrahedron (for center of mass estimation), notice that
        % the last one is calculated as the estimated total mass minus the
        % rest mass parameters
        ctrllnputVal(5:8)=[unCertainParam1(8:10); unCertainParam1(1)-sum(unCertainParam1(8:10))];
        
        % Get System Properties
        [~,~,~,~,~,Quat,MM,~,GFF,GFU,JacU,JacCons1,~,~,CenMat,InertLeftComponent,InertRightComponent1]...
            =System_TAWE_mex(tSpan1(ii),simStateVal,ctrlParamVal,ctrllnputVal,nhSignal);
        MM=sum(MM,3); % Sum up inertia matrix of all bodies
        JacU=sum(JacU,3).'; % Sum up constraint Jacobian Matrices of all inputs
        JacCons1=sum(JacCons1,3); % Sum up constraint Jacobian Matrices of all constraints
        GFF=sum(GFF,2); % Sum up generalized forces (excluding those from inputs and load tetrahedron)
        GFU=sum(GFU(:,1:4),2); % Sum up load tetrahedron generalized forces (equivalent to estimated generalized uncertainload)
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
            System_TAWE_mex(tSpan1(ii),simStateVal2,ctrlParamVal,ctrllnputVal,nhSignal2);
        
        
        
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
        loadRotMat=quat2Mat(Quat(:,5)); % Frame rotation of the uncertain body
        
        % Intermeditate terms
        midTerm1=refddt-KK1Gain*(simStateVal(10:11)-refdt); % derivative of zeta in the paper (Eq.(15))
        midTerm2=refdt-KK1Gain*(simStateVal(1:2)-ref); % zeta in the paper (Eq.(15))
        midTerm3=simStateVal(10:11) - refdt + KK1Gain*(simStateVal(1:2)-ref); % xi in the paper (Eq.(15))
        
        % Control Input using estimated model properties
        JJww=0.5*JacU2Tpose.'; % estimated disturbance jacobian matrix
        RRInv=(JJww.'*JJww+KK2Gain); % R from the paper (Eq.(15))
        
        
        % The controller design from Eq.(14-16) using estimated model parameters
        taweInput=(JacU2Tpose)\( MMCombine*midTerm1... % compensate inertia acceleration 
                       + CCCombine*midTerm2 ... % compensate coriolis force
                       - GFFCombine... % compensate generalized force
                       - JlamMapFull.'*GFU... % compensate uncertain load (from load tetrahedron (for center of mass estimation))
                       + JlamMapFull.'*(unCertainParam2.*qqdot)... % compensate uncertain damping
                       - c1Coeff * RRInv * midTerm3); % feedback term (Eq.(15))
        
        %Prepare for Uncertain Parameter Jacobian Calculation
        uncertainInertRightComponent=InertRightComponent(1:6,:)*JlamMapFull; % Right inertia component matrix 
        uncertainInertRightComponentDot=InertRightComponentDot(1:6,:)*JlamMapFull...
                             +InertRightComponent(1:6,:)*JlamMapFullDot; % Right inertia component matrix derivative
        
        % Jacobian for parameter from uncertain body inertia and load
        paramJac1 = uncertainBodyParamJacobianFunc_TAWE_mex...
                    (simStateVal(10:11),midTerm1,midTerm2,loadRotMat,...
                    uncertainInertRightComponent,uncertainInertRightComponentDot,(JacU3Tpose.'));
        % Jacobian for parameter from uncertain damping
        paramJac2=diag(qqdot)*JlamMapFull;
                
        % Update uncertain parameters, note that (simFreq/ctrlFreq) 
        % ensures the step size is the same for different control sample 
        % rates
        unCertainParam1 = unCertainParam1 ...
                        - hh * GammaInv1 * paramJac1 * midTerm3 * (simFreq/ctrlFreq); 
        unCertainParam2 = unCertainParam2 ...
                        - hh * GammaInv2 * paramJac2 * midTerm3 * (simFreq/ctrlFreq);
        
    end
    
    simInputVal(3:4)=taweInput; % Insert into input for simulation
    
    % Using Runge Kutta 4th method to integrate system ODE
    [simStateVal,nhSignal,consForce]=rk4(@Flow_TAWE_mex,hh,flowNum,tSpan1(ii),...
                                        simStateVal,simParamVal,simInputVal,nhSignal,consForce);
    
    % The codes below refreshes the 3D visualization for run-time animation
    % Uncomment to run.
%     if rem(ii,15)==0 %adjust the frame rate through the number
%         Sim.drawNow(tSpan1(ii),simStateVal,simParamVal,simInputVal,nhSignal); % requires 3D environment
%     end
end

toc;
save(strcat('./Data/simTAWEDataNoWKI'),'simStateData'); % Save result

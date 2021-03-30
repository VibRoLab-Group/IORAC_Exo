%%  No-Disturbance Stationary Upper Limb Exoskeleton (Arm4D) simulations
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script contains the code for stationary upper limb exoskeleton (Arm4D) simulations
% under no disturbance. It contains the simulations with two uncertainty conditions 
% (i.e., involving Link 4 uncertainty only, and involving both Link 3 and Link 4 uncertainties).


PathSetup; % Include folders to path directory and reset workspace

%% 3D Visualization (Uncomment to run)
% Below is the initialization codes for the 3D visualization of the
% stationary exoskeleton, which includes 3D models, load tetrahedrons, and
% main coordinate frames. Note that the human body models are not included.
% Uncomment the code in this section to initiate visualization. (Requires
% OpenGL for smooth animation)

%{
% Create figure
Sim=sAxes('Arm 4-DOF',3,'Frames_Arm4D.mat',@numTF_Arm4D_mex);
Sim.setAxesProp('BaseFrame',[-0.3 -0.1 -0.35; 0.3 1.2 0.3],[135 30]).setPresetTF('WORLD');

% Load coordinate frame markers
[MainLinkChain,COMChain,link3PTChain,link4PTChain]...
    =Sim.genPlot({'MainLinkChain';'COMChain';'link3PTChain';'link4PTChain'});
COMChain.setLineSpec('none','*','b',5).setPlotPoint({'Link1COMFrame';'Link2COMFrame';'Link3COMFrame';'Link4COMFrame'});
link3PTChain.setLineSpec(':','+','g',3).setPlotPoint({'Link3mpt1Frame';'Link3mpt2Frame';'Link3mpt3Frame';'Link3mpt4Frame';'Link3mpt2Frame';'Link3mpt1Frame';'Link3mpt3Frame';'Link3mpt1Frame';'Link3mpt4Frame';});
link4PTChain.setLineSpec('--','x','r',3).setPlotPoint({'Link4mpt1Frame';'Link4mpt2Frame';'Link4mpt3Frame';'Link4mpt4Frame';'Link4mpt2Frame';'Link4mpt1Frame';'Link4mpt3Frame';'Link4mpt1Frame';'Link4mpt4Frame';});

% Load .stl mesh files
[BasePatch,Link1Patch,Link2Patch,Link3Patch,Link4Patch]...
=Sim.genPatch({'BasePatch' 'BaseFrame';...
               'Link1Patch' 'Link1Frame';...
               'Link2Patch' 'Link2Frame';...
               'Link3Patch' 'Link3Frame';...
               'Link4Patch' 'Link4Frame';...
               });
BasePatch.setFaceProp([0.8 0.8 1],0.75).setModel('arm4DBasePart.STL',1,[],0.0005);
Link1Patch.setFaceProp('c',0.75).setModel('arm4DLink1.STL',1,[],0.0005);
Link2Patch.setFaceProp('b',0.75).setModel('arm4DLink2.STL',1,[],0.0005);
Link3Patch.setFaceProp('g',0.75).setModel('arm4DLink3.STL',1,[],0.0005);
Link4Patch.setFaceProp('r',0.75).setModel('arm4DLink4.STL',1,[],0.0005);
grid on;
%}

%% Select uncertainty condition
% The number determines the uncertainty case (Link 4 uncertainty only, or Link 3 and Link 4 uncertainties)
iiTest=2; % Set to 1 for Link 4 uncertainty, set to 2 for Link 3 and 4 uncertainties

%% Load Model Parameters
TestNumList=[4 34]; 
Arm4D_TestNum=TestNumList(iiTest); % see above for options

Arm4D_Parameters; % Load model parameters (depends on Arm4D_TestNum)
paramSelectVec=1:14; % Uncertain parameter ranges
if Arm4D_TestNum==4
    paramSelectVec=[5:8,12:14]; % For Link 4 uncertain parameters only
elseif Arm4D_TestNum==34
    tSpan = 0:hh:600; % Longer simulation time span for this case
end

load(strcat('simLink',num2str(Arm4D_TestNum),'Ref'),'refData','refdtData','refddtData');

tic;

%% Simulation Variables


% simulation state (first half are joint angles, and second half are angular velocities)
simStateVal=[initialPositionVal;zeros(4,1)]; 
simParamVal=ppSimVal; % Simulation uses true parameters

% The uncertain model paramters, which will be updated
ppCtrlVal2Update=[ppInertiaUnknownVal;ppLoadUnknwonVal];
% Controller uses estimated parameters
ctrlParamVal=[ppCtrlVal1;ppCtrlVal2Update(1:8);ppCtrlVal3]; 


% The input vectors uu are defined as 
% (1-4) Torque input at the exoskeleton joints
% (5-8) Link 3 virtual input at the load tetrahedron vertices (-z direction) 
%       (for center of mass estimation), used to acquire the Jacobian for 
%       uncertain load parameter more conveniently. These inputs are not 
%       used for simulation, but for controller design only.
% (9-12) Link 4 virtual input at the load tetrahedron vertices (-z direction) 
%       (for center of mass estimation), used to acquire the Jacobian for 
%       uncertain load parameter more conveniently. These inputs are not 
%       used for simulation, but for controller design only.
simInputVal=zeros(12,1); % System input for simulation
ctrlInputVal=zeros(12,1); % This is used for controller design only 

jointInput=zeros(4,1); % Control input at the four exoskeleton joints

consForce=0; % Constraint force (none for this case)
nhSignal=0; % Nonholonomic constraint (none for this case)

%% Simulation Loop
simStateData=zeros(8,numel(tSpan)); % Stores trajectory tracking performance
simParamData=zeros(14,numel(tSpan)); % Stores parameters estimation performance

for ii=1:numel(tSpan)
    
    tt=tSpan(ii); % Current time variable
    simStateData(:,ii)=simStateVal; 
    simParamData(:,ii)=ppCtrlVal2Update;
    
    % Current tracking reference
    ref=refData(:,ii);
    refdt=refdtData(:,ii);
    refddt=refddtData(:,ii);
    
    if rem(ii,simFreq/ctrlFreq)==0 % used for control sample rate
        
        % Get System Properties
        [~,~,~,~,~,Quat,MM,~,GFF,GFU,JacU,~,~,~,CenMat,InertLeftComponent,InertRightComponent]=...
                System_Arm4D_mex(tSpan(ii),simStateVal,ctrlParamVal,ctrlInputVal,nhSignal);
        MM=sum(MM,3); % Sum up inertia matrix of all bodies
        JacU=sum(JacU,3).'; % Sum up constraint Jacobian Matrices of all inputs
        GFF=sum(GFF,2); % Sum up generalized forces (excluding those from inputs and load tetrahedron)
        CenMat=sum(CenMat,3); % Sum up all centripetal matrices'
        
        % Sum up load tetrahedron generalized forces (equivalent to estimated 
        % generalized uncertain loads). For the case that only Link 4 has
        % uncertainty, since the link 3 load tetrahedron parameters are
        % zeros, only loads from link 4 load tetrahedron exist.
        GFU=sum(GFU(:,1:8),2); 
        
        % Properties used to numerically calculate Coriolis matrix and 
        % constraint Jacobian derivative
        simStateValNextStep=simStateVal;
        % Use the current state velocity to estimate the position after a
        % short time step
        simStateValNextStep(1:4)=simStateValNextStep(1:4)+simStateVal(5:8)*(hh);
        % Obtain estimated properties after a short time step
        [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,InertRightComponentNextStep]=...
            System_Arm4D_mex(tSpan(ii),simStateValNextStep,ctrlParamVal,ctrlInputVal,nhSignal);
        
        % Numerically calculate coriolis matrix
        InertRightComponentDot=(InertRightComponentNextStep-InertRightComponent)/(hh);
        CorMat=InertLeftComponent*InertRightComponentDot;
        
        % Intermeditate terms
        midTerm1=refddt-KK1Gain*(simStateVal(5:8)-refdt); % derivative of zeta in the paper (Eq.(15))
        midTerm2=refdt-KK1Gain*(simStateVal(1:4)-ref); % zeta in the paper (Eq.(15))
        midTerm3=simStateVal(5:8) - refdt + KK1Gain*(simStateVal(1:4)-ref); % xi in the paper (Eq.(15))
        
        % Control input and parameter update
        JacUInput=JacU(1:4,:); % Jacobian of exoskeleton joint control input 
        JacULink3=JacU(5:8,:); % Jacobian of Link 3 load tetrahedron forces
        JacULink4=JacU(9:12,:); % Jacobian of Link 4 load tetrahedron forces
        rotMatLink3=quat2Mat(Quat(:,5)); % Rotation of Link 3 Frame
        rotMatLink4=quat2Mat(Quat(:,6)); % Rotation of Link 4 Frame
        
        c1Coeff=4; % c1Coeff is overwrited to 4 for the no-disturbance simulations
        JJww=JacUInput; % Estimated disturbance jacobian matrix
        RRInv=(JJww.'*JJww+KK2Gain); % R from the paper (Eq.(15))
        jointInput=(JacUInput.')\( MM*midTerm1...% compensate inertia acceleration 
                       + (CorMat + CenMat)*midTerm2 ... % compensate coriolis force
                       - GFF... % compensate generalized force
                       - GFU - ... % compensate uncertain load (from load tetrahedron (for center of mass estimation))
                       c1Coeff * RRInv * midTerm3); % feedback term (Eq.(15))

                   
        
        %Prepare for Link 4 Uncertain Parameter Jacobian Calculation
        uncertainInertRightComponentLink4=InertRightComponent(31:36,:);  % Right inertia component matrix for Link 4 uncertainty
        uncertainInertRightComponentDotLink4=InertRightComponentDot(31:36,:); % Derivative of Right inertia component matrix for Link 4 uncertainty
        
        
        % Jacobian for parameter from Link 4 uncertain body inertia and load
        paramJacLink4 = uncertainBodyParamJacobianFunc_Arm4D_mex...
                        (simStateVal(5:8),midTerm1,midTerm2,rotMatLink4,uncertainInertRightComponentLink4,...
                        uncertainInertRightComponentDotLink4,JacULink4);
        if Arm4D_TestNum==34
            %Prepare for Link 3 Uncertain Parameter Jacobian Calculation
            uncertainInertRightComponentLink3=InertRightComponent(25:30,:); % Right inertia component matrix for Link 4 uncertainty
            uncertainInertRightComponentDotLink3=InertRightComponentDot(25:30,:); % Derivative of Right inertia component matrix for Link 4 uncertainty
            
            % Jacobian for parameter from Link 3 uncertain body inertia and load
            paramJacLink3 = uncertainBodyParamJacobianFunc_Arm4D_mex...
                        (simStateVal(5:8),midTerm1,midTerm2,rotMatLink3,uncertainInertRightComponentLink3,...
                        uncertainInertRightComponentDotLink3,JacULink3);
                    
            % Parameter Jacobian for all uncertain parameters
            paramJac=[paramJacLink3(1:4,:);paramJacLink4(1:4,:);paramJacLink3(5:7,:);paramJacLink4(5:7,:);];
        else
            paramJac=paramJacLink4; % For link 4 uncertainty only
        end
        
        % Update selected parameters based on uncertainty condition
        ppCtrlVal2Update(paramSelectVec) = ppCtrlVal2Update(paramSelectVec)...
                                           - hh * GammaInv * paramJac * midTerm3 * (simFreq/ctrlFreq);
                                
        % Insert updated parameters for upcoming controller design                     
        if Arm4D_TestNum==34
            % Insert estimated mass parameters on Link 3 load tetrahedron 
            % (for center of mass estimation), notice that the last one
            % is calculated the estimated total mass minus the rest mass 
            % parameters
            ctrlInputVal(5:8)=[ppCtrlVal2Update(9:11); ppCtrlVal2Update(1)-sum(ppCtrlVal2Update(9:11))];
        end
        % Insert estimated mass parameters on Link 4 load tetrahedron 
        % (for center of mass estimation), notice that the last one
        % is calculated as the estimated total mass minus the rest mass 
        % parameters
        ctrlInputVal(9:12)=[ppCtrlVal2Update(12:14); ppCtrlVal2Update(5)-sum(ppCtrlVal2Update(12:14))];
        ctrlParamVal=[ppCtrlVal1;ppCtrlVal2Update(1:8);ppCtrlVal3];
    end
    
    
    simInputVal(1:4)=jointInput; % Insert joint control input to simulation
    % Using Runge Kutta 4th method to integrate system ODE
    [simStateVal,nhSignal,consForce]=rk4(@Flow_Arm4D_mex,hh,flowNum,tSpan(ii),...
                                    simStateVal,simParamVal,simInputVal,nhSignal,consForce);

    % The codes below refreshes the 3D visualization for run-time animation
    % Uncomment to run.
%     if rem(ii,15)==0 %adjust the frame rate through the number
%         Sim.drawNow(tSpan(ii),simStateVal,simParamVal,simInputVal,nhSignal); % requires 3D environment
%     end
end

toc;
save(strcat('./Data/simLink',num2str(Arm4D_TestNum),'Data'),'simStateData','simParamData'); % Save result

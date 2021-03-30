%%  Stationary Upper Limb Exoskeleton (Arm4D) simulations under Disturbance
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script contains the code for Stationary Upper Limb (Arm4D) exoskeleton simulations
% under disturbance. It contains both Link 3 and Link 4 uncertainties. The
% performance of PD, SMC, and IO-RAC feedback controllers are generated for
% comparison


PathSetup; % Include folders to path directory and reset workspace

%% Note
% In the paper, three controllers are tuned with respect to the sample
% reference trajectory, so that they yield similar performance level. 
% However, since references and disturbances are randomly generated, the 
% current tuned controllers may yield different performance levels for
% another set of data. However, the qualitative observation about each
% controller remains the same.


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

%% Load Model Parameters
        
Arm4D_TestNum=34; % Both link 3 and link 4 has uncertainties
Arm4D_Parameters; % Load model parameters (depends on Arm4D_TestNum)
paramSelectVec=1:14; % Uncertain parameter ranges

tSpan=0:hh:240; % Simulation time span

load('simLink34RefDisturbance','refData','refdtData','refddtData','disturbData');

%% Simulation Loop

% These matrix stores results
simStateData=zeros(8,numel(tSpan),3); % for tracking performance
inputData=zeros(4,numel(tSpan),3); % for feedback input

% kk=1 for sliding mode controller (SMC)
% kk=2 for PD controller (High Gain)
% kk=3 for proposed IO-RAC
tic;
for kk=1:3


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
    feedbackVal=zeros(4,1); % feedback control input
    gainSwitchVal=zeros(4,1); % gain switching control input
    
    consForce=0; % Constraint force (none for this case)
    nhSignal=0; % Nonholonomic constraint (none for this case)
    
    % Stores result for the current loop
    simStateDataSingle=zeros(8,numel(tSpan)); 
    inputDataSingle=zeros(4,numel(tSpan));
    
    for ii=1:numel(tSpan)
        
        tt=tSpan(ii);
        simStateDataSingle(:,ii)=simStateVal;
        ref=refData(:,ii);
        refdt=refdtData(:,ii);
        refddt=refddtData(:,ii);

        if rem(ii,simFreq/ctrlFreq)==0
            
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

            JJww=0.4*JacUInput; % Estimated disturbance jacobian matrix
            JJww=JJww*diag([3;1;3;1]);
            if kk==1
                RRInv=(0*JJww.'*JJww+KK2Gain);
            elseif kk==2
                RRInv=(0*JJww.'*JJww+4*KK2Gain);
            else
                RRInv=(JJww.'*JJww+KK2Gain);
            end
            feedbackVal = - (JacUInput.')\(c1Coeff * RRInv * midTerm3); % feedback term (Eq.(15))
            jointInput=    feedbackVal + ...
                           (JacUInput.')\( MM*midTerm1...% compensate inertia acceleration 
                           + (CorMat + CenMat)*midTerm2 ... % compensate coriolis force
                           - GFF... % compensate generalized force
                           - GFU); % compensate uncertain load (from load tetrahedron (for center of mass estimation))
            if kk==1
                % The following gain switching contoller uses the true sign function
                % gainSwitchVal = -(JacU.')\(diag(0.5*[2;1;2;1])*sign(midTerm3)); 

                % The following gain switching contoller uses the mimic
                % sign function 2*(sigmoid(scaleMid3*x)-0.5), which is
                % designed for slightly reduce chattering
                scaleMid3 = 1e2*midTerm3; % a large scaling factor 
                gainSwitchVal = - (JacUInput.')\(diag(0.5*[2;1;2;1])*(1-exp(-scaleMid3))./(1+exp(-scaleMid3)));
                jointInput = jointInput + gainSwitchVal;
            end
            
                            
            %Prepare for Link 4 Uncertain Parameter Jacobian Calculation
            uncertainInertRightComponentLink4=InertRightComponent(31:36,:);  % Right inertia component matrix for Link 4 uncertainty
            uncertainInertRightComponentDotLink4=InertRightComponentDot(31:36,:); % Derivative of Right inertia component matrix for Link 4 uncertainty
            % Jacobian for parameter from Link 4 uncertain body inertia and load
            paramJacLink4 = uncertainBodyParamJacobianFunc_Arm4D_mex...
                            (simStateVal(5:8),midTerm1,midTerm2,rotMatLink4,uncertainInertRightComponentLink4,...
                            uncertainInertRightComponentDotLink4,JacULink4);
            %Prepare for Link 3 Uncertain Parameter Jacobian Calculation
            uncertainInertRightComponentLink3=InertRightComponent(25:30,:); % Right inertia component matrix for Link 4 uncertainty
            uncertainInertRightComponentDotLink3=InertRightComponentDot(25:30,:); % Derivative of Right inertia component matrix for Link 4 uncertainty
            % Jacobian for parameter from Link 3 uncertain body inertia and load
            paramJacLink3 = uncertainBodyParamJacobianFunc_Arm4D_mex...
                        (simStateVal(5:8),midTerm1,midTerm2,rotMatLink3,uncertainInertRightComponentLink3,...
                        uncertainInertRightComponentDotLink3,JacULink3);

            % Parameter Jacobian for all uncertain parameters
            paramJac=[paramJacLink3(1:4,:);paramJacLink4(1:4,:);paramJacLink3(5:7,:);paramJacLink4(5:7,:);];
            % Update selected parameters based on uncertainty condition
            ppCtrlVal2Update(paramSelectVec) = ppCtrlVal2Update(paramSelectVec)...
                                               - hh * GammaInv * paramJac * midTerm3 * (simFreq/ctrlFreq);
  
            % Insert estimated mass parameters on Link 3 load tetrahedron 
            % (for center of mass estimation), notice that the last one
            % is calculated the estimated total mass minus the rest mass 
            % parameters
            ctrlInputVal(5:8)=[ppCtrlVal2Update(9:11); ppCtrlVal2Update(1)-sum(ppCtrlVal2Update(9:11))];
            % Insert estimated mass parameters on Link 4 load tetrahedron 
            % (for center of mass estimation), notice that the last one
            % is calculated as the estimated total mass minus the rest mass 
            % parameters
            ctrlInputVal(9:12)=[ppCtrlVal2Update(12:14); ppCtrlVal2Update(5)-sum(ppCtrlVal2Update(12:14))];
            ctrlParamVal=[ppCtrlVal1;ppCtrlVal2Update(1:8);ppCtrlVal3];
        end
        
        inputDataSingle(:,ii) = feedbackVal; % stores the current feedback term
        if kk==1
            inputDataSingle(:,ii) = feedbackVal + gainSwitchVal; % include the gain switching term for the SMC case
        end
        
        simInputVal(1:4) = jointInput + disturbData(:,ii).*[3;1;3;1]; % Insert joint input overlaid with noise
        % Using Runge Kutta 4th method to integrate system ODE
        [simStateVal,nhSignal,consForce]=rk4(@Flow_Arm4D_mex,hh,flowNum,tSpan(ii),simStateVal,simParamVal,simInputVal,nhSignal,consForce);
        
        % The codes below refreshes the 3D visualization for run-time animation
        % Uncomment to run.
    %     if rem(ii,15)==0 %adjust the frame rate through the number
    %         Sim.drawNow(tSpan(ii),simStateVal,simParamVal,simInputVal,nhSignal); % requires 3D environment
    %     end
    end
    
    % collect data from all cases
    simStateData(:,:,kk)=simStateDataSingle; 
    inputData(:,:,kk)=inputDataSingle;
end

toc;
save(strcat('./Data/simLink',num2str(Arm4D_TestNum),'DataDisturbance'),'simStateData','inputData');
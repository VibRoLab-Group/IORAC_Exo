%% Uncertain Parameter Jacobian Function Generation for the Wearable Wrist Exoskeleton (TAWE) simulations
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script generates the uncertain parameter Jacobian function
% introduced by a unknown body. The uncertain parameters includes mass,
% moment of inertia, and load parameters at load tetrahedron (for center 
% of mass estimation). The generated function is used in wearable wrist
% exoskeleton (TAWE) simulations.

PathSetup;
PathSetup; % Reset path twice to clean cache

%% Symbolic Formulation of Jacobian Function

% Uncertain local inertia parameters: (1) Mass; (2-7) Moments of Inertia (xx,xy,xz,yy,yz,zz)
uncertainSym=sym('unknownInert%d',[7,1]);

% Rotation of the local frame
rotMat=sym('rotMat%d',[3,3]);
% Inertia matrix in global frame
MM=[eye(3)*uncertainSym(1) zeros(3);
    zeros(3) rotMat*formInertiaMatrix_(uncertainSym(2:end))*rotMat.'];

vv1=sym('v1Anony%d',[2,1]); % term that represents the derivative of zeta (Eq.(15))
vv2=sym('v2Anony%d',[2,1]); % term that represents zeta
qqdot=sym('qqdot%d',[2,1]); % velocity of the state (of the constrained model) 

% Right inertia component matrix
inertRightComponentMat=sym('inertRightComponentMatSym%d',[6,2],'real');
% The derivative of right inertia component matrix
inertRightComponentMatDot=sym('inertRightComponentMatDotSym%d',[6,2],'real'); 
wJac=inertRightComponentMat(4:6,:); % the part respect to angular velocity
angVel=wJac*qqdot; % angular velocity
% left inertia component matrix
inertLeftComponentMat=inertRightComponentMat.'*MM;

% The acceleration related to inertia matrix
inertTerm1=inertLeftComponentMat*inertRightComponentMat*vv1;
% The acceleration related to coriolis matrix
corTerm1=inertLeftComponentMat*inertRightComponentMatDot*vv2;
% The acceleration related to centripetal matrix
cenTerm1=wJac.'*skew3(angVel)*MM(4:6,4:6)*wJac*vv2;
% The total acceleration
allTerm1=inertTerm1+corTerm1+cenTerm1;

% Symbolically calculate jacobian of uncertainSym from the total acceleration
paramJac1=jacobian(allTerm1,uncertainSym).';
% The input jacobian for mass parameter on the load tetrahedron (for center 
% of mass estimation)
paramJac2=sym('paramJac2%d',[4,2],'real');

% Since sum of mass on load tetrahedron equals to the total mass, a
% transformation matrix is designed as
transformMat=[eye(10); 1 zeros(1,6) -ones(1,3)];
% Combine the two jacobians (Notice the negative sign is because it is on 
% the right side of the equation in the simulation multibody model)
paramJac=[paramJac1;-paramJac2];
% Transform the matix yields the final jacobian
paramJac=transformMat.'*paramJac; 

%% Function Generation
% Generate function script
funcName='./ModelFunction_TAWE/uncertainBodyParamJacobianFunc_TAWE';
paramJacFunc=matlabFunction(paramJac,'file',funcName,'vars',{qqdot,vv1,vv2,...
            rotMat,inertRightComponentMat,inertRightComponentMatDot,paramJac2});

% Generate precompiled functions (comment this part if MATLAB compiler is not availble)
codegen(funcName,'-d','./ModelFunction_TAWE/','-args',{double(qqdot*0),double(vv1*0),...
        double(vv2*0),double(rotMat*0),double(inertRightComponentMat*0),...
        double(inertRightComponentMatDot*0),double(paramJac2*0)},'-jit');
delete("./*.mexw64");
delete("./*.mexa64");


%% Local Function
% form inertia matrix (angular part) from parameters
function outVar_=formInertiaMatrix_(inVar_)
    diagTerm_=inVar_([1 4 6]);
    offDiagMat_=sym(0)*zeros(3,3);
    offDiagMat_(1,2)=inVar_(2);
    offDiagMat_(1,3)=inVar_(3);
    offDiagMat_(2,3)=inVar_(5);
    
    outVar_ = diag(diagTerm_) + offDiagMat_ + offDiagMat_.';
end
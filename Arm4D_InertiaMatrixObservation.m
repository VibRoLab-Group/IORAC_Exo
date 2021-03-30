%% Observation of the Inertia Matrix of the Stationary Upper Limb (Arm4D) Exoskeleton
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script demonstrate how Eq.(22) is obtained and shows that it is not  
% uncertain inertia parameters, but their linear combinations that uniquely
% affect the inertia matrix.

PathSetup; % Include folders to path directory
load ModelInfo_Arm4D.mat; % Load the symbolic information of the model

%% Symbolic Operation

% This part eliminates the known inertial properties. The remaining part is
% the model uncertainty.
inertialMatrix = ModelInfo.Dynamics.InertialMatrix;
knownVariables = ModelInfo.Variables.Discrete([1:28,56:58]); 
% Note that the Ang0 variables are base orientation angles of the 
% exoskeleton, whose values are zeros (i.e., no default orientation)
unknownInertialMatrix = subs(inertialMatrix,knownVariables,zeros(31,1)); 

% Observe that there are eight uncertain inertia variables/parameters
uncertainInertialVariable = symvar(unknownInertialMatrix);
uncertainInertialVariable = uncertainInertialVariable(10:17).';
momentTerms = uncertainInertialVariable([1:3,5:7]); %The unknown moment of inertia terms

% Create Symbolic Variable for Eq.(22)
momentLinearComboSym = sym('momentLinearCombo',[3,1]);
momentLinearComboVal = [
                        momentTerms(5) - momentTerms(6);
                        momentTerms(2) + momentTerms(6);
                        momentTerms(3) + momentTerms(6);
                      ]; % Expression of Eq.(22)

% Here, we will substitute some of the unknown moments with the expression
% in Eq.(22)
momentSubSym = [momentTerms(5); momentTerms(2); momentTerms(3)]; % Terms to be substituted
momentSubVal = momentLinearComboSym + momentSubSym - momentLinearComboVal; % Expressions of the substituted terms
unknownInertialMatrixSub=subs(unknownInertialMatrix,momentSubSym,momentSubVal);
unknownInertialMatrixSub=simplify(expand(unknownInertialMatrixSub)); % Symbolic simplification

% After the substitution, notice that there are only five remaining moment 
% of inertia terms. The term "linkUkn4i33" has been canceled out in the 
% inertial matrix.
remainingMomentTerm=intersect([momentTerms;momentLinearComboSym],symvar(unknownInertialMatrixSub));
disp(remainingMomentTerm.');

% Hence, there are redundancy in the uncertain inertial parameters. Even
% through these parameters exist in the multibody model, the inertial
% matrix is uniquely affected by "remainingMomentTerm", but not every
% variables in "momentTerms".

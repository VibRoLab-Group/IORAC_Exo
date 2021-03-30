function [modelState]=ekfStep(modelFunc,modelState)
% An iteration step of extended Kalman filter process.
% Details of the process are commented below.

% (No Copyright Claimed)
% Author: Jiamin Wang; Revised: 2-Aug-2020

    pp=modelState.pp; % Parameters of the model (will not be updated)
    xx=modelState.xx; % State of the model (posteriori of the last step)
    uu=modelState.uu; % Input of the model (e.g., sensor measurement, control input)
    PPxx=modelState.PPxx; % Covariance matrix of the state (posteriori of the last step)
    QQ=modelState.QQ; % Covariance of process noise
    RR=modelState.RR; % Covariance of observation noise
    
    QQBase=modelState.QQBase; % Auxiliary linear covariance of process noise to avoid singularity
    RRBase=modelState.RRBase; % Auxiliary linear covariance of observation noise to avoid singularity

    xxHat=modelFunc.xxPlus(pp,xx,uu); % Priori state estimation
    Fxx=modelFunc.Fxx(pp,xx,uu); % Calculate the state Jacobian from the state model df/dx
    Fvv=modelFunc.Fvv(pp,xx,uu); % Calculate the noise Jacobian from the state model df/dv
    
    modelState.yyHat=modelFunc.yy(pp,xxHat,uu); % observation based on priori State
    Hxx=modelFunc.Hxx(pp,xxHat,uu); % Calculate the state Jacobian from the observation model df/dx
    Hww=modelFunc.Hww(pp,xxHat,uu); % Calculate the noise Jacobian from the observation model df/dw
    
    PPxxHat=Fxx*PPxx*Fxx.'+Fvv*QQ*Fvv.'+QQBase; % Priori State covariance Estimation
    PPyy=Hxx*PPxxHat*Hxx.'+Hww*RR*Hww.'+RRBase; % Observation covariance Estimation
    PPxy=PPxxHat*Hxx.'; % covariance between state and observation
    
    GG=PPxy/PPyy; % Optimal Kalman gain
    xx=xxHat-GG*modelState.yyHat; % Posteriori state update
    modelState.xx=modelFunc.xxPost(pp,xx,uu); % Post-processing to eliminate numerical errors
    modelState.PPxx=PPxxHat-GG*PPyy*GG.'; % Posteriori state covariance update
    
end
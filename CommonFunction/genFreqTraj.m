function [trajOut0, trajOut1, trajOut2]=genFreqTraj(inDim,inTime,inAmp,inFreq)
% Generate periodic or quasiperiodic trajectories (and their 1st and 2nd
% order time derivatives) composed by harmonic waves of different frequencies, 
% amplitudes, and phases within a provided time span. Note that the
% amplitudes are randomly generated.
% Outputs:
% (1) trajOut0: the generated trajectories
% (2) trajOut1: the generated trajectories (1st order time derivative)
% (3) trajOut2: the generated trajectories (2nd order time derivative)
% Inputs:
% (1) inOrder: the dimensionality of the trajectories (rows)
% (2) inTime: the time span (columns)
% (3) inAmp: the base magnitude of the amplitudes
% (4) inFreq: the frequency components of the trajectories

% (No Copyright Claimed)
% Author: Jiamin Wang; Revised: 2-Aug-2020

freqNum=numel(inFreq); 
tSpan=reshape(inTime,1,[]);
dataNum=numel(tSpan);
freqArray=2*reshape(inFreq,[],1)*pi;

trajOut0=zeros(inDim,dataNum); 
trajOut1=zeros(inDim,dataNum);
trajOut2=zeros(inDim,dataNum);

    for ii=1:inDim
        
        ampArray=inAmp.*(rand([freqNum,2])-0.5); %This randomly generates the amplitudes for the sine and cosine waves
        
        traj0=(ampArray(:,1).')*([sin(freqArray.*tSpan)])...
            +(ampArray(:,2).')*([cos(freqArray.*tSpan)]);
        traj1=(ampArray(:,1).')*(freqArray.*[cos(freqArray.*tSpan)])...
            -(ampArray(:,2).')*(freqArray.*[sin(freqArray.*tSpan)]);
        traj2=-(ampArray(:,1).')*(freqArray.*freqArray.*[sin(freqArray.*tSpan)])...
            -(ampArray(:,2).')*(freqArray.*freqArray.*[cos(freqArray.*tSpan)]);
        
        trajOut0(ii,:)=traj0;
        trajOut1(ii,:)=traj1;
        trajOut2(ii,:)=traj2;

    end

end
%% Plot Generations for the Wearable Wrist Exoskeleton (TAWE) simulations
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script generates the plots for the wearable wrist exoskeleton 
% (TAWE) simulations. The plots generated are for Fig. 6 in the paper.

PathSetup; % Include folders to path directory and reset workspace
drawingParam; % Drawing parameters for plot generations

%% Generate Fig. 6 in the paper
% This section plots the figure for the wearable wrist Exoskeleton (TAWE) 
% simulations.

TAWE_Parameters; % Load TAWE parameters
tSpan=0:hh:60;

% Load references and simulated trajectories
load('simTAWERef.mat');
load(strcat('simTAWEDataNoWKI'),'simStateData');

% Generate plots
Fig=figure(435);
set(Fig, 'Units', 'pixels', 'OuterPosition', [0, 0, 576*2, 514*2]);
delete(Fig.Children)

jj=1;
axes;
set(gca,'Position',subplotPosSet(3,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan,refData(1,1:numel(tSpan)),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold on
plot(tSpan,simStateData(1,1:numel(tSpan)),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 60])
ylim([-1.1 1.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\theta_\mathrm{RUD}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('$$r_{\theta,\mathrm{RUD}}$$','$$\theta_\mathrm{RUD}$$','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=2;
axes;
set(gca,'Position',subplotPosSet(3,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan,refData(2,1:numel(tSpan)),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold on
plot(tSpan,simStateData(2,1:numel(tSpan)),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 60])
ylim([-1.1 1.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\theta_\mathrm{FE}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('$$r_{\theta,\mathrm{FE}}$$','$$\theta_\mathrm{FE}$$','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=3;
axes;
set(gca,'Position',subplotPosSet(3,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan,simStateData(1,1:numel(tSpan))-refData(1,1:numel(tSpan)),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold on
plot(tSpan,simStateData(2,1:numel(tSpan))-refData(2,1:numel(tSpan)),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 60])
ylim([-0.05 0.05])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\mathbf{\epsilon}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('$$\epsilon_{\theta,\mathrm{RUD}}$$','$$\epsilon_{\theta,\mathrm{FE}}$$','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

load(strcat('simTAWEDataWKI'),'simStateData');

tSpan=0:hh:120;

for jj=4:6
    
axes;
set(gca,'Position',subplotPosSet(3,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan,simStateData(1,1:numel(tSpan),jj-3)-refData(1,1:numel(tSpan)),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold on
plot(tSpan,simStateData(2,1:numel(tSpan),jj-3)-refData(2,1:numel(tSpan)),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
xlim([0 120])
ylim([-0.1 0.1])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\mathbf{\epsilon}$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('$$\epsilon_{\theta,\mathrm{RUD}}$$','$$\epsilon_{\theta,\mathrm{FE}}$$','interpreter','latex','FontSize',26,'orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

end
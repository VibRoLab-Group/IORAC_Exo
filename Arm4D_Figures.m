%% Plot Generations for the Stationary Upper Limb Exoskeleton simulations
% Simulation code for "Inverse Optimal Robust Adaptive Controller for 
% Upper Limb Rehabilitation Exoskeletons with Inertia and Load Uncertainties"
% by Jiamin Wang (jmechw@vt.edu) and Oumar R. Barry (obarry@vt.edu)
% Deparment of Mechanical Engineering, Virginia Tech
% 
% This script generates the plots for the stationary upper limb exoskeleton
% simulations. The plots generated are for Fig. 2, Fig. 3, and Fig. 4 in
% the paper.

PathSetup;
drawingParam;

%% Generate Fig. 2 in the paper
% This section plots the figure for the stationary upper limb exoskeleton
% simulation case, which assumes no disturbance, and only Link 4 has model
% uncertainties.

TestNumList=[4 34];
iiTest=1;
Arm4D_TestNum=TestNumList(iiTest);
Arm4D_Parameters;
paramSelectVec=1:14;
if Arm4D_TestNum==4
    paramSelectVec=[5:8,12:14];
end
tSpan=0:hh:300;

% Load references and simulated trajectories
load(strcat('simLink',num2str(Arm4D_TestNum),'Ref'),'refData','refdtData','refddtData');
load(strcat('simLink',num2str(Arm4D_TestNum),'Data'));

% Generate plots
Fig=figure(431);
set(Fig, 'Units', 'pixels', 'OuterPosition', [0, 0, 576*2, 514*2.5]);
delete(Fig.Children)

for jj=1:4
    axes;
    set(gca,'Position',subplotPosSet(4,2,jj,0.3,0.15,0.05,0.05));
    plot(tSpan,simStateData(jj,1:numel(tSpan))-initialPositionVal(jj),'-r','linewidth',3);
    hold on
    plot(tSpan,refData(jj,1:numel(tSpan))-initialPositionVal(jj),'--k','linewidth',2);
    xlim([0 60])
    hold off; grid on;
    xlabel('\textbf{Time (s)}','interpreter','latex')
    ylabel(strcat('\textbf{$$\theta_',num2str(jj),'$$  (rad)}'),'interpreter','latex')
    set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
    legend('\textbf{traj.}','\textbf{ref.}','interpreter','latex','orientation','horizontal','location','south')
    text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
    grid on;
end

jj=5;
axes;
set(gca,'Position',subplotPosSet(4,2,jj,0.3,0.15,0.05,0.05));
for kk=1:4
    if kk==2
        hold on
    end
    plot(tSpan,simStateData(kk,1:numel(tSpan)) - refData(kk,1:numel(tSpan)),lineStyleDefault{kk},'linewidth',lineWidthDefault(kk),'Color',lineColorDefault{kk});
end
xlim([0 60])
ylim([-0.4 0.5])
hold off; grid on;
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel(strcat('\textbf{$$\mathbf{\epsilon}_\theta$$  (rad)}'),'interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('\textbf{$$\theta_1$$}','\textbf{$$\theta_2$$}','\textbf{$$\theta_3$$}','\textbf{$$\theta_4$$}','interpreter','latex','orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=6;
axes;
set(gca,'Position',subplotPosSet(4,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan,simParamData(5,:)-CheckVal1(5),'b','linewidth',3)
hold on
ttspan1=[0 60];
hold off; grid on;
xlim(ttspan1)
ylim([-8 2])
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel('\textbf{$$\tilde{m}_2$$ (kg)}','interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

jj=7;
axes;
set(gca,'Position',subplotPosSet(4,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan,simParamData(6,:)-CheckVal1(6),lineStyleDefault{1},'linewidth',lineWidthDefault(1),'Color',lineColorDefault{1})
hold on
ttspan1=[0 300];
plot(tSpan,simParamData(7,:)-CheckVal1(7),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
plot(tSpan,simParamData(8,:)-CheckVal1(8),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold off; grid on;
xlim(ttspan1)
ylim([-1 0.5])
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel('\textbf{$$\tilde{\mathbf{\phi}}_2$$ (kg$$\cdot$$m$$^2$$)}','interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('\textbf{$${\phi}_{2,xx}$$}','\textbf{$${\phi}_{2,yy}$$}','\textbf{$${\phi}_{2,zz}$$}','interpreter','latex','orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

link4UncertainCOMEst = (tretraVertices4 * [simParamData(12:end,:);simParamData(5,:)-sum(simParamData(12:end,:),1)])./simParamData(5,:);

jj=8;
axes;
set(gca,'Position',subplotPosSet(4,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan,link4UncertainCOMEst(1,:)-link4UncertainCOM(1),lineStyleDefault{1},'linewidth',lineWidthDefault(1),'Color',lineColorDefault{1})
hold on
ttspan1=[0 300];
plot(tSpan,link4UncertainCOMEst(2,:)-link4UncertainCOM(2),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
plot(tSpan,link4UncertainCOMEst(3,:)-link4UncertainCOM(3),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold off; grid on;
xlim(ttspan1)
ylim([-0.1 0.2])
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel('\textbf{$$\tilde{\mathbf{d}}_2$$ (m)}','interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('\textbf{$${\mathbf{d}_{2,x}}$$}','\textbf{$${\mathbf{d}_{2,y}}$$}','\textbf{$${\mathbf{d}_{2,z}}$$}','interpreter','latex','orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;


%% Generate Fig. 3 in the paper
% This section plots the figure for the stationary upper limb exoskeleton
% simulation case, which assumes no disturbance, and both Link 3 and Link 4 
% has model uncertainties.

TestNumList=[4 34];
iiTest=2;
Arm4D_TestNum=TestNumList(iiTest);
Arm4D_Parameters;
paramSelectVec=1:14;
if Arm4D_TestNum==4
    paramSelectVec=[5:8,12:14];
end

tSpan=0:hh:600;

% Load references and simulated trajectories
load(strcat('simLink',num2str(Arm4D_TestNum),'Ref'),'refData','refdtData','refddtData');
load(strcat('simLink',num2str(Arm4D_TestNum),'Data'));

% Generate plots
Fig=figure(432);
set(Fig, 'Units', 'pixels', 'OuterPosition', [0, 0, 576*2, 514*1.4]);
delete(Fig.Children)

jj=1;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan,simParamData(1,:)-CheckVal1(1),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
hold on
ttspan1=[0 600];
plot(tSpan,simParamData(5,:)-CheckVal1(5),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold off; grid on;
xlim(ttspan1)
ylim([-2 2])
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel('\textbf{$$\tilde{m}$$ (kg)}','interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('\textbf{$$m_1$$}','\textbf{$$m_2$$}','interpreter','latex','orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;


jj=2;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan,simParamData(1,:)-CheckVal1(1),lineStyleDefault{1},'linewidth',lineWidthDefault(1),'Color',lineColorDefault{1})
hold on
ttspan1=[0 600];
plot(tSpan,simParamData(5,:)-CheckVal1(5),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
plot(tSpan,simParamData(8,:)-CheckVal1(8),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold off; grid on;
xlim(ttspan1)
ylim([-2 2])
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel('\textbf{$$\tilde{\mathbf{\phi}}$$ (kg$$\cdot$$m$$^2$$)}','interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('\textbf{$${\phi}_{1,xx}$$}','\textbf{$${\phi}_{2,xx}$$}','\textbf{$${\phi}_{2,zz}$$}','interpreter','latex','orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;

linCombVal=CheckVal1([7 3 4])+[-CheckVal1(8);CheckVal1(8);CheckVal1(8)];
jj=3;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
plot(tSpan,simParamData(7,:)-simParamData(8,:)-linCombVal(1),lineStyleDefault{1},'linewidth',lineWidthDefault(1),'Color',lineColorDefault{1})
hold on
ttspan1=[0 600];
plot(tSpan,simParamData(3,:)+simParamData(8,:)-linCombVal(2),lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
plot(tSpan,simParamData(4,:)+simParamData(8,:)-linCombVal(3),lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold off; grid on;
xlim(ttspan1)
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel('\textbf{$$\tilde{\mathbf{c}}_\phi$$ (kg$$\cdot$$m$$^2$$)}','interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('\textbf{$${c}_{\phi,1}$$}','\textbf{$${c}_{\phi,2}$$}','\textbf{$${c}_{\phi,3}$$}','interpreter','latex','orientation','horizontal','location','south')
grid on;
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');


link3UncertainCOMEst = (tretraVertices3 * [simParamData(9:11,:);simParamData(1,:)-sum(simParamData(9:11,:),1)])./simParamData(1,:);
link4UncertainCOMEst = (tretraVertices4 * [simParamData(12:end,:);simParamData(5,:)-sum(simParamData(12:end,:),1)])./simParamData(5,:);
link3UncertainCOMErrNorm = sum((link3UncertainCOMEst-link3UncertainCOM).^2,1).^0.5;
link4UncertainCOMErrNorm = sum((link4UncertainCOMEst-link4UncertainCOM).^2,1).^0.5;

jj=4;
axes;
set(gca,'Position',subplotPosSet(2,2,jj,0.3,0.15,0.05,0.05));
ttspan1=[0 600];
plot(tSpan,link3UncertainCOMErrNorm,lineStyleDefault{2},'linewidth',lineWidthDefault(2),'Color',lineColorDefault{2})
hold on
plot(tSpan,link4UncertainCOMErrNorm,lineStyleDefault{3},'linewidth',lineWidthDefault(3),'Color',lineColorDefault{3})
hold off; grid on;
xlim(ttspan1)
ylim([-0.1 0.2])
xlabel('\textbf{Time (s)}','interpreter','latex')
ylabel('\textbf{$$\| \tilde{\mathbf{d}} \|$$ (m)}','interpreter','latex')
set(gca,'fontname','times new roman','FontSize',26,'FontWeight','bold')
legend('\textbf{$${\mathbf{d}_{1}}$$}','\textbf{$${\mathbf{d}_{2}}$$}','interpreter','latex','orientation','horizontal','location','south')
text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',26,'Color','black','fontname','times new roman');
grid on;


%% Generate Fig. 4 in the paper
% This section plots the figure for the stationary upper limb exoskeleton
% simulation case, which assumes model disturbance, and both Link 3 and Link 4 
% has model uncertainties. The performances of PD, SMC, and proposed IO-RAC
% are compared in this figure.

TestNumList=[4 34];
iiTest=2;
Arm4D_TestNum=TestNumList(iiTest);
Arm4D_Parameters;
paramSelectVec=1:14;
if Arm4D_TestNum==4
    paramSelectVec=[5:8,12:14];
end
tSpan=0:hh:240;

% Load references and simulated trajectories
load(strcat('simLink',num2str(Arm4D_TestNum),'RefDisturbance'),'refData','refdtData','refddtData');
load(strcat('simLink',num2str(Arm4D_TestNum),'DataDisturbance'));

% Generate plots
Fig=figure(433);
set(Fig, 'Units', 'pixels', 'OuterPosition', [0, 0, 576*2, 514*2.5]);
delete(Fig.Children)

tTicks=1:numel(tSpan);
for jj=1:4
    axes;    
    set(gca,'Position',subplotPosSet(4,2,jj,0.3,0.15,0.05,0.05));
    hold on;
    if rem(jj,2)==1
        plt2=plot(tSpan,simStateData(jj,:,2) - refData(jj,tTicks),lineStyleDefault{2},'linewidth',3.5,'Color',lineColorDefault{2});
        plt3=plot(tSpan,simStateData(jj,:,3) - refData(jj,tTicks),lineStyleDefault{3},'linewidth',2.5,'Color',lineColorDefault{3});
    else
        plt3=plot(tSpan,simStateData(jj,:,3) - refData(jj,tTicks),lineStyleDefault{3},'linewidth',2.5,'Color',lineColorDefault{3});
        plt2=plot(tSpan,simStateData(jj,:,2) - refData(jj,tTicks),lineStyleDefault{2},'linewidth',3.5,'Color',lineColorDefault{2});
    end
    plt1=plot(tSpan,simStateData(jj,:,1) - refData(jj,tTicks),lineStyleDefault{1},'linewidth',2.5,'Color',lineColorDefault{1});
    hold off; grid on;
    xlim([120 150])
    ylim([-0.125 0.125])
    if jj==1
        ylim([-0.1 0.1])
    end
    xlabel('\textbf{Time (s)}','interpreter','latex')
    ylabel(strcat('\textbf{$$\epsilon_{\theta,',num2str(jj),'}$$  (rad)}'),'interpreter','latex')
    legend([plt2,plt1,plt3],{'\textbf{PD}','\textbf{SMC}','\textbf{IO-RAC}'},'interpreter','latex','orientation','horizontal','location','south')
    set(gca,'fontname','times new roman','FontSize',24,'FontWeight','bold')
    text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',24,'Color','black','fontname','times new roman');
    grid on;
end

tTicks=1:numel(tSpan);
for jj=5:8
    axes;
    set(gca,'Position',subplotPosSet(4,2,jj,0.3,0.15,0.05,0.05));
    hold on;
    if jj~=7
        plt1=plot(tSpan,inputData(jj-4,:,1) ,lineStyleDefault{1},'linewidth',2,'Color',lineColorDefault{1});
    end
    if rem(jj,2)==0
        plt2=plot(tSpan,inputData(jj-4,:,2) ,lineStyleDefault{2},'linewidth',2,'Color',lineColorDefault{2});
        plt3=plot(tSpan,inputData(jj-4,:,3) ,lineStyleDefault{3},'linewidth',2,'Color',lineColorDefault{3});
        hold off; grid on;
        ylim([-2 1])
    else
        plt3=plot(tSpan,inputData(jj-4,:,3) ,lineStyleDefault{3},'linewidth',2,'Color',lineColorDefault{3});
        if jj==7
            plt1=plot(tSpan,inputData(jj-4,:,1) ,'--','linewidth',2,'Color',lineColorDefault{1});
        end
        plt2=plot(tSpan,inputData(jj-4,:,2) ,lineStyleDefault{2},'linewidth',2,'Color',lineColorDefault{2});
        hold off; grid on;
        ylim([-4 2])
        if jj==5
            ylim([-3 2])
        end
    end
    xlim([120 150])
    xlabel('\textbf{Time (s)}','interpreter','latex')
    ylabel(strcat('\textbf{$$u_{b,\theta,',num2str(jj-4),'}$$  (N$$\cdot$$m)}'),'interpreter','latex')
    legend([plt2,plt1,plt3],{'\textbf{PD}','\textbf{SMC}','\textbf{IO-RAC}'},'interpreter','latex','orientation','horizontal','location','south')
    set(gca,'fontname','times new roman','FontSize',24,'FontWeight','bold')
    text(get(gca,'xlim')*[0.075;0.925],get(gca,'ylim')*[0.1;0.9],lineLabelDefault{jj},'horizontalalignment','center','FontWeight','bold','FontSize',24,'Color','black','fontname','times new roman');
    grid on;
end

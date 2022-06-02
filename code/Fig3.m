%% Figures @NREM_comp
%
% Theoretical and Computational Neuroscience (TCN), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 3. Numerical integration pipeline
%
% Created: Mon 1 Feb 2021, 10:46
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Sun 22 May 2022, 17:19
% Last edited by: Martin Esparza-Iaizzo
%

%% Definitions

% r = mean firing rate
% a = adaptation

% Arbitrary units for the time constants in order to avoid dimensions
tau_r = 1;
tau_a = 25; 

% Fixed values for sigmoidal function
x0 = 5; 
r0 = 0.5;
k = 15;

% Initial conditions
ii_cc = [0;0];

% Noise function is employed with Ornstein - Uhlenbeck method

duration = 400; %s
save_dt = 0.01; % Precision on saving values in noise function and in solving the ode
dt = 0.01; % must be the same as save_dt
theta = 0.05; % Noise "amplitude"
sigma = 0.25; % Sigma1
numsignals = 1;

% Generate noise realization. 
[noise, noise_t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals,6);

% Select regime
regime = 1; % 1 = Bistable or 0 = Oscillatory'

if regime == 1
    b = 1;
    w = 6.3;
    I = 2.35;
end
if regime == 0
    b = 1;
    w = 6;
    I = 2.5;
end

% Integrate model. 
[t,ra] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0), [0:dt:duration], ii_cc);
r1 = ra(:,1)';


regime = 0; % 1 = Bistable or 0 = Oscillatory'
ii_cc = [0.5;0.5];
if regime == 1
    b = 1;
    w = 6.3;
    I = 2.35;
end
if regime == 0
    b = 1;
    w = 6;
    I = 2.5;
end

[t1,ra1] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0), [0:dt:duration], ii_cc);
r1 = ra1(:,1)';

%% Generate noise windows. 
start = 50; window_duration = duration - start;
new_noise = noise(1:start/dt);
complete_window(:,1) = noise;
for i = 2:10
    [window,window_t] = OUNoiseWindow(theta,0.5,window_duration,dt,save_dt,numsignals,noise(start/dt+1),29);
    complete_window(:,i) = [new_noise; window];
end


%% Plotting 

f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 15 17]*1.75);

fontSize = 16;

% Sigmas.
sigma1 = 0.25; % Sigma 1
[noise1, noise_t1] = OUNoise(theta,sigma1,duration,dt,save_dt,numsignals,8); % Generate noise realization

ax = axes('Position',[0.095 0.86 0.25 0.075]);
ax.PositionConstraint = 'innerposition';
plot(t,noise1,'LineWidth',.5,'Color',[0.5 0.5 0.5]);
ylim([-1 1])
xticks([]); xticklabels({''}); 
xlabel('Time (AU)')
set(ax,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it \xi(t)}','FontName','Arial')

sigma2 = 0.02; % Sigma 2
[noise1, noise_t1] = OUNoise(theta,sigma2,duration,dt,save_dt,numsignals,8); % Generate noise realization

ax = axes('Position',[0.475 0.86 0.25 0.075]);
ax.PositionConstraint = 'innerposition';
plot(t,noise1,'LineWidth',.5,'Color',[0.5 0.5 0.5]);
ylim([-0.1 0.1])
xticks([]); xticklabels({''}); 
ylabel('{\it \xi(t)}')
xlabel('Time (AU)')
set(ax,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it \xi(t)}','FontName','Arial')

% Scenarios
% Up Bistable
ax1_1 = axes('Position',[0.075 0.45 0.175 0.15]);
ax1_1.PositionConstraint = 'innerposition';
plot(t,r1,'LineWidth',2.5,'Color','#027EDC'); hold on; 
xlim([0 175])
ylim([-1 1.1])
xticks([20 50]); xticklabels({'0','10'}); 
yticks([]); yticklabels({''});
xlabel('Time (AU)')
alpha(ax1_1,0)
set(ax1_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it            r(t)}','FontName','Arial')

ax1_2 = axes('Position',[0.075 0.45 0.171 0.06]);
ax1_2.PositionConstraint = 'innerposition';
plot(t,noise,'LineWidth',.5,'Color','k');
xlim([0 175])
ylim([-0.5 0.75])
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
alpha(0)
set(ax1_2,'FontSize',fontSize,'Box','off','LineWidth',1.5,'FontName','Arial')
ylabel('{\it\xi(t)}','FontName','Arial')

% Down Bistable
ax1_1 = axes('Position',[0.075 0.62 0.175 0.15]);
ax1_1.PositionConstraint = 'innerposition';
plot(t,r1,'LineWidth',2.5,'Color','#027EDC'); hold on; 
xlim([100 300])
ylim([-1 1.1])
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
set(ax1_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it            r(t)}','FontName','Arial')

ax1_2 = axes('Position',[0.075 0.62 0.171 0.06]);
ax1_2.PositionConstraint = 'innerposition';
plot(t,noise,'LineWidth',.5,'Color','k');
xlim([100 300])
ylim([-0.75 0.5])
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
set(ax1_2,'FontSize',fontSize,'Box','off','LineWidth',1.5,'FontName','Arial')
ylabel('{\it\xi(t)}','FontName','Arial')

% Up Oscillatory
ax1_1 = axes('Position',[0.075 0.22 0.175 0.15]);
ax1_1.PositionConstraint = 'innerposition';
plot(t1,r1,'LineWidth',2.5,'Color','#FF44C8'); hold on; 
xlim([200 375])
ylim([-1 1.1])
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
set(ax1_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it            r(t)}','FontName','Arial')

ax1_2 = axes('Position',[0.075 0.22 0.171 0.06]);
ax1_2.PositionConstraint = 'innerposition';
plot(t1,noise,'LineWidth',.5,'Color','k');
xlim([200 375])
ylim([-0.75 0.5])
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
set(ax1_2,'FontSize',fontSize,'Box','off','LineWidth',1.5,'FontName','Arial')
ylabel('{\it\xi(t)}','FontName','Arial')

% Down Oscillatory
ax1_1 = axes('Position',[0.075 0.05 0.175 0.15]);
ax1_1.PositionConstraint = 'innerposition';
plot(t1,r1,'LineWidth',2.5,'Color','#FF44C8'); hold on; 
xlim([50 300])
ylim([-1 1.1])
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
set(ax1_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it            r(t)}','FontName','Arial')

ax1_2 = axes('Position',[0.075 0.05 0.171 0.06]);
ax1_2.PositionConstraint = 'innerposition';
plot(t1,noise,'LineWidth',.5,'Color','k');
xlim([50 300])
ylim([-0.75 0.5])
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
xlabel('Time (AU)')
set(ax1_2,'FontSize',fontSize,'Box','off','LineWidth',1.5,'FontName','Arial')
ylabel('{\it\xi(t)}','FontName','Arial')


% Noise
distances = [0.075 0.175 0.275 0.49 0.59 0.69];
distances = fliplr(distances);
realizations = [1 2 3 298 299 300];
for i = 1:5
    ax2_2 = axes('Position',[0.375 distances(i) 0.225 0.08]);
    ax2_2.PositionConstraint = 'innerposition';
    plot(t,complete_window(:,i+1),'LineWidth',.6,'Color','k'); hold on; 
    xticks([]); xticklabels({''}); 
    yticks([-1 1]); yticklabels({'-1','1'});
    ylim([-1.25 1.25])
    str = sprintf('_{%i}',realizations(i));
    str = ['{\it\xi(t)}',str];
    xlim([0 175])
    set(ax2_2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')
    ylabel(str,'FontName','Arial')
end
ax2_3 = axes('Position',[0.375 distances(6) 0.225 0.08]);
ax2_3.PositionConstraint = 'innerposition';
plot(t,complete_window(:,i+1),'LineWidth',.6,'Color','k'); hold on; 
xticks([20 50]); xticklabels({'0','10'}); 
yticks([-1 1]); yticklabels({'-1','1'});
ylim([-1.25 1.25])
xlabel('Time (AU)')
xlim([0 175])
set(ax2_3,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('\xi(t)_{300}','FontName','Arial')

% Stimuli
% Stimuli 1
ax3_1 = axes('Position',[0.725 0.69 0.23 0.08]);
ax3_1.PositionConstraint = 'innerposition';
plot(t,[zeros(1,5000) -1.25*ones(1,3000) zeros(1,32001)],'LineWidth',2.5,'Color','#DAA520');
xticks([]); xticklabels({''}); 
yticks([-0.75 0.75]); yticklabels({'-50','50'});
ylim([-1.5 1.75])
xlim([0 175])
set(ax3_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

ax3_2 = axes('Position',[0.865 0.735 0.08 0.03]);
ax3_2.PositionConstraint = 'innerposition';
plot(t,complete_window(:,3),'LineWidth',0.25,'Color','k');
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
ylim([-1.5 1.5])
xlim([0 175])
set(ax3_2,'FontSize',fontSize-2,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it \xi(t)}_2','FontName','Arial'); 

% Stimuli 2
ax3_1 = axes('Position',[0.725 0.59 0.23 0.08]);
ax3_1.PositionConstraint = 'innerposition';
plot(t,[zeros(1,5000) -1*ones(1,3000) zeros(1,32001)],'LineWidth',2.5,'Color','#DAA520');
xticks([]); xticklabels({''}); 
yticks([-0.75 0.75]); yticklabels({'-50','50'});
ylim([-1.5 1.75])
xlim([0 175])
set(ax3_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

ax3_2 = axes('Position',[0.865 0.635 0.08 0.03]);
ax3_2.PositionConstraint = 'innerposition';
plot(t,complete_window(:,3),'LineWidth',0.25,'Color','k');
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
ylim([-1.5 1.5])
xlim([0 175])
set(ax3_2,'FontSize',fontSize-2,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it \xi(t)}_2','FontName','Arial'); 

% Stimuli 3
ax3_1 = axes('Position',[0.725 0.49 0.23 0.08]);
ax3_1.PositionConstraint = 'innerposition';
plot(t,[zeros(1,5000) -0.75*ones(1,3000) zeros(1,32001)],'LineWidth',2.5,'Color','#DAA520');
xticks([]); xticklabels({}); 
yticks([-0.75 0.75]); yticklabels({'-50','50'});
ylim([-1.5 1.75])
xlim([0 175])
set(ax3_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

ax3_2 = axes('Position',[0.865 0.535 0.08 0.03]);
ax3_2.PositionConstraint = 'innerposition';
plot(t,complete_window(:,3),'LineWidth',0.25,'Color','k');
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
ylim([-1.5 1.5])
xlim([0 175])
set(ax3_2,'FontSize',fontSize-2,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it \xi(t)}_2','FontName','Arial'); 

% Stimuli 4
ax3_1 = axes('Position',[0.725 0.275 0.23 0.08]);
ax3_1.PositionConstraint = 'innerposition';
plot(t,[zeros(1,5000) 0.75*ones(1,3000) zeros(1,32001)],'LineWidth',2.5,'Color','#DAA520');
xticks([]); xticklabels({''}); 
yticks([-0.75 0.75]); yticklabels({'-50','50'});
% ylabel('Hz')
ylim([-1.5 1.75])
xlim([0 175])
set(ax3_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

ax3_2 = axes('Position',[0.865 0.32 0.08 0.03]);
ax3_2.PositionConstraint = 'innerposition';
plot(t,complete_window(:,3),'LineWidth',0.25,'Color','k');
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
ylim([-1.5 1.5])
xlim([0 175])
set(ax3_2,'FontSize',fontSize-2,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it \xi(t)}_2','FontName','Arial'); 

% Stimuli 5
ax3_1 = axes('Position',[0.725 0.175 0.23 0.08]);
ax3_1.PositionConstraint = 'innerposition';
plot(t,[zeros(1,5000) 1*ones(1,3000) zeros(1,32001)],'LineWidth',2.5,'Color','#DAA520');
xticks([]); xticklabels({''}); 
yticks([-0.75 0.75]); yticklabels({'-50','50'});
% ylabel('Hz')
ylim([-1.5 1.75])
xlim([0 175])
set(ax3_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

ax3_2 = axes('Position',[0.865 0.22 0.08 0.03]);
ax3_2.PositionConstraint = 'innerposition';
plot(t,complete_window(:,3),'LineWidth',0.25,'Color','k');
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
ylim([-1.5 1.5])
xlim([0 175])
set(ax3_2,'FontSize',fontSize-2,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it\xi(t)}_2','FontName','Arial'); 

% Stimuli 6
ax3_1 = axes('Position',[0.725 0.075 0.23 0.08]);
ax3_1.PositionConstraint = 'innerposition';
plot(t,[zeros(1,5000) 1.25*ones(1,3000) zeros(1,32001)],'LineWidth',2.5,'Color','#DAA520');
xticks([20 50]); xticklabels({'0','10'}); 
yticks([-0.75 0.75]); yticklabels({'-50','50'});
 xlabel('Time (AU)')
ylim([-1.5 1.75])
xlim([0 175])
set(ax3_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

ax3_2 = axes('Position',[0.865 0.12 0.08 0.03]);
ax3_2.PositionConstraint = 'innerposition';
plot(t,complete_window(:,3),'LineWidth',0.25,'Color','k');
xticks([]); xticklabels({''}); 
yticks([]); yticklabels({''});
ylim([-1.5 1.5])
xlim([0 175])
set(ax3_2,'FontSize',fontSize-2,'Box','on','LineWidth',1.5,'FontName','Arial')
ylabel('{\it \xi(t)}_2','FontName','Arial'); 

%% Save figure ––––––– Uncomment and edit to save to personalised location

% cd '/Volumes/GoogleDrive-101271366273470520077/My Drive/PaperBelen/Figures/Temp figures'
% set(f,'Renderer','Painter')
% exportgraphics(gcf,'fig3.pdf','Resolution',300,'BackgroundColor','none')
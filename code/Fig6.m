%% Figures @NREM_comp
%
% Theoretical and Computational Neuroscience (TCN), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 6. Firing rate response to changes in stimulus firing frequency
%
% Created: Mon 8 Mar 2021, 12:31
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Mon 31 May 2022, 13:30
% Last edited by: Martin Esparza-Iaizzo

%% Definitions.

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
ii_cc_up = [0;0];
ii_cc_down = [1;1];
% Noise function is employed with Ornstein - Uhlenbeck method

duration = 750; %s
save_dt = 0.01; % Precision on saving values in noise function and in solving the ode
dt = 0.01; % must be the same as save_dt
theta =0.05; %About comparable...
sigma = 0.25;
numsignals = 1;

[noise_up, noise_t_up] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals,2);
[noise_down, noise_t_down] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals,17);
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

%% GENERATE SPIKE TRAIN

% Time stamps for UPs and DOWNs
start_1_up = 325; 
start_2_up = 200;

start_1_down = 575;
start_2_down = 400;

duration_s = 40; % stimuli duration
freqs = [10:10:50]; % Incremental frequencies
spike_train = generateStimuli(duration_s,dt,1,freqs);

stimuli_up = zeros(1,duration/dt+1); % Up stimuli
stimuli_0 = zeros(1,duration/dt+1); % No stimuli
stimuli_up(start_1_up/dt:(start_1_up+duration_s)/dt-1) = spike_train(10,:);
stimuli_up(start_2_up/dt:(start_2_up+duration_s)/dt-1) = spike_train(1,:);

stimuli_down = zeros(1,duration/dt+1); % Down stimuli
stimuli_down(start_1_down/dt:(start_1_down+duration_s)/dt-1) = spike_train(10,:);
stimuli_down(start_2_down/dt:(start_2_down+duration_s)/dt-1) = spike_train(1,:);

%% Simulation

% Integrate the model
[t,ra] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise_up,noise_t_up,x0,k,r0), [0:dt:duration], ii_cc_up);
[t,ra_s1] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,noise_up,noise_t_up,x0,k,r0,stimuli_up), [0:dt:duration], ii_cc_up);
[t,ra_s2_0] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise_down,noise_t_down,x0,k,r0), [0:dt:duration], ii_cc_down);
[t,ra_s2] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,noise_down,noise_t_down,x0,k,r0,stimuli_down), [0:dt:duration], ii_cc_down);

r = ra(:,1)';
r_s1 = ra_s1(:,1)';
r_s2_0 = ra_s2_0(:,1)';
r_s2 = ra_s2(:,1)';

%% Load example data

addpath(genpath('/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/variables/sigma = 0.25'))
load('post_up_bis.mat')
load('post_up_osc.mat')
load('post_down_bis.mat')
load('post_down_osc.mat')


%% Redefine parameters for second part of the figure

% Arbitrary units for the time constants in order to avoid dimensions
tau_r = 1;
tau_a = 25; 

% Fixed values for sigmoidal function
x0 = 5; 
r0 = 0.5;
k = 15;

% Initial conditions
ii_cc_up = [0;0];

% Noise function is employed with Ornstein - Uhlenbeck method

duration = 400; %s
save_dt = 0.01; % Precision on saving values in noise function and in solving the ode
dt = 0.01; % must be the same as save_dt
theta =0.05; %About comparable...
sigma = 0.25;
numsignals = 1;

[noise2, noise_t2] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals,6);

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

%% Simulation

% Integrate the model
[t2,ra2] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise2,noise_t2,x0,k,r0), [0:dt:duration], ii_cc_up);
r2 = ra2(:,1)';

%% Plotting

f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17 20]*1.75);
fontSize = 16;
LineW = 2.5; 

% A1
ax1 = axes('Position',[0.11 0.85 0.375 0.12]);
ax1.PositionConstraint = 'innerposition';
plot(t,r,'LineWidth',LineW,'Color','k'); hold on; 
plot(t,stimuli_0,'LineWidth',LineW,'Color','#DAA520')
ylim([-1 1])
xlim([100 500])
xticklabels({''})
yticks([-1 0 1])
yticklabels({'-1 ','0 ','1 '})
ylabel('{\it r(t) / \color[rgb]{.85 .647 .125}s(t)}')
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')


% A2
ax2 = axes('Position',[0.11 0.7 0.375 0.12]);
ax2.PositionConstraint = 'innerposition';
plot(t,r_s1,'LineWidth',LineW,'Color','k'); hold on; 
plot(t,stimuli_up*0.9,'LineWidth',LineW,'Color','#DAA520')
ylim([-1 1])
xlim([100 500])
yticks([-1 0 1])
yticklabels({'-1 ','0 ','1 '})
xlabel('Time (AU)')
ylabel('{\it r(t) / \color[rgb]{.85 .647 .125}s(t)}')
set(ax2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

% A3
ax3 = axes('Position',[0.575 0.7 0.35 0.27]);
ax3.PositionConstraint = 'innerposition';
plot([-100:10:100],post_up_osc(12,:),'.','MarkerSize',45,'Color','#027EDC','LineWidth',1)
ylim([0 20])
yticks([0 5 10 15 20])
yticklabels({'0','5','10','15','20'})
ylabel('Duration (AU)')
xlabel('Frequency (Hz)')
set(ax3,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

% B1
ax1 = axes('Position',[0.11 0.475 0.375 0.12]);
ax1.PositionConstraint = 'innerposition';
plot(t,r_s2_0,'LineWidth',LineW,'Color','k'); hold on; 
plot(t,stimuli_0,'LineWidth',LineW,'Color','#DAA520')
ylim([-1 1])
xlim([300 700])
xticklabels({''})
yticks([-1 0 1])
yticklabels({'-1 ','0 ','1 '})
ylabel('{\it r(t) / \color[rgb]{.85 .647 .125}s(t)}')
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')


% B2
ax2 = axes('Position',[0.11 0.325 0.375 0.12]);
ax2.PositionConstraint = 'innerposition';
plot(t,r_s2,'LineWidth',LineW,'Color','k'); hold on; 
plot(t,stimuli_down*0.9,'LineWidth',LineW,'Color','#DAA520')
ylim([-1 1])
xlim([300 700])
yticks([-1 0 1])
yticklabels({'-1 ','0 ','1 '})
xlabel('Time (AU)')
ylabel('{\it r(t) / \color[rgb]{.85 .647 .125}s(t)}')
set(ax2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')


% B3
ax4 = axes('Position',[0.575 0.325 0.35 0.27]);
ax4.PositionConstraint = 'innerposition';
plot([-100:10:100],post_down_osc(15,:),'.','MarkerSize',45,'Color','#027EDC','LineWidth',1)
ylim([0 20])
yticks([0 5 10 15 20])
yticklabels({'0','5','10','15','20'})
ylabel('Duration (AU)')
xlabel('Frequency (Hz)')
set(ax4,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')


positions_inset_x = 0.025;
positions_inset_y = 0.025;

ax5 = axes('Position',[0.575+positions_inset_x 0.325+positions_inset_y 0.125 0.1]);
ax5.PositionConstraint = 'innerposition';
plot(t2,r2,'LineWidth',1.75,'Color','k');
xlim([100 300])
ylim([-0.1 1.1])
xticks([110 140])
xticklabels({'0','10'})
yticklabels({''})
yticks([])
set(ax5,'FontSize',fontSize-2,'Box','on','LineWidth',1.5,'FontName','Arial')

positions_inset_x = 0.025;
positions_inset_y = 0.15;

ax6 = axes('Position',[0.575+positions_inset_x 0.7+positions_inset_y 0.125 0.1]);
ax6.PositionConstraint = 'innerposition';
plot(t2,r2,'LineWidth',1.75,'Color','k');
xlim([0 200])
ylim([-0.1 1.1])
xticks([0 35])
xticklabels({'0','10'})
yticks([])
yticklabels({''})
set(ax6,'FontSize',fontSize-2,'Box','on','LineWidth',1.5,'FontName','Arial')

%% Save figure ––––––– Uncomment and edit to save to personalised location

% cd '/Volumes/GoogleDrive-101271366273470520077/My Drive/PaperBelen/Figures/Temp figures'
% set(f,'Renderer','Painter')
% exportgraphics(gcf,'fig6.pdf','Resolution',300,'BackgroundColor','none')





%% Figures @NREM_comp
%
% Theoretical and Computational Neuroscience (TCN), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 4. Description of the interaction between noise and firing rate
%
% Created: Sat 27 Feb 2021, 11:57
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Mon 23 May 2022, 12:04
% Last edited by: Martin Esparza-Iaizzo
%
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

duration = 1500; %s
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


%% Simulation

% Integrate model
[t,ra] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0), [0:dt:duration], ii_cc);
r_A = ra(:,1)';
a_A = ra(:,2)';

start = 50; % Start of simulation
new_noise = noise(1:start/dt); window_duration = duration - start;

[window_B,window_t_B] = OUNoiseWindow(theta,sigma,window_duration,dt,save_dt,numsignals,noise(start/dt+1),1);
complete_window_B = [new_noise; window_B];

[window_C,window_t_C] = OUNoiseWindow(theta,sigma,window_duration,dt,save_dt,numsignals,noise(start/dt+1),22);
complete_window_C = [new_noise; window_C];

[t,ra_B] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,complete_window_B,noise_t,x0,k,r0), [0:dt:duration], ii_cc);
[t,ra_C] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,complete_window_C,noise_t,x0,k,r0), [0:dt:duration], ii_cc);

r_B = ra_B(:,1)';
r_C = ra_C(:,1)';
%% Plotting FIG1
f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17 20]*1.75);

fontSize = 16;

ax = axes('Position',[0.1 0.8 0.85 0.15]);
ax.PositionConstraint = 'innerposition';
plot(t,noise,'LineWidth',1.5,'Color','r'); hold on 
plot(t,r_A,'LineWidth',2.5,'Color','k'); 
ylim([-1.1 1.1])
yticks([-1 0 1]); yticklabels({'-1','0','1'})
xlabel('Time (AU)')
ylabel('{\it r(t) / \color{red}\xi(t)}')
set(ax,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')


ax1 = axes('Position',[0.1 0.375 0.85 0.15]);
ax1.PositionConstraint = 'innerposition';
plot(t,complete_window_B,'LineWidth',1.5,'Color','r'); hold on;
plot(t,complete_window_C,'LineWidth',1.5,'Color',[1 .7 .7]);
ylim([-1.1 1.1])
xlim([0 450])
yticks([-1 0 1]); yticklabels({'-1','0','1'})
xticks([0 150 300 450 600]); xticklabels({'0','150','300','450','600'})
xlabel('Time (AU)')
ylabel('{\it \xi(t)}')
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')


ax2 = axes('Position',[0.1 0.55 0.85 0.15]);
ax2.PositionConstraint = 'innerposition';
plot(t,r_B,'LineWidth',3.25,'Color','k'); hold on;
plot(t,r_C,'LineWidth',3.25,'Color',[.7 .7 .7]);
ylim([-0.1 1.1])
xlim([0 450])
yticks([-1 0 1]); yticklabels({'-1','0','1'})
xticks([0 150 300 450 600]); xticklabels({''})
ylabel('{\it r(t)}')
set(ax2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

%% Save figure ––––––– Uncomment and edit to save to personalised location

% cd '/Volumes/GoogleDrive-101271366273470520077/My Drive/PaperBelen/Figures/Temp figures'
% set(f,'Renderer','Painter')
% exportgraphics(gcf,'fig4.pdf','Resolution',300,'BackgroundColor','none')

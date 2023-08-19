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
addpath(genpath('./'))

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
ii_cc_down = [1;0];
% Noise function is employed with Ornstein - Uhlenbeck method

duration = 750; %s
save_dt = 0.01; % Precision on saving values in noise function and in solving the ode
dt = 0.01; % must be the same as save_dt
theta =0.05; %About comparable...
sigma = 0.25;
numsignals = 1;

[noise_up, noise_t_up] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals,2);
[noise_down, noise_t_down] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals,100);
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
start_up = 175;

start_down = 92;

duration_s = 10; % stimuli duration
freqs = [10:10:100]; % Incremental frequencies
spike_train = generateStimuli(duration_s,dt,1,freqs);

stimuli = zeros(1,duration/dt+1); % Stimuli

% Noise params up
window_duration_up = duration - start_up;
new_noise_up = noise_up(1:start_up/dt);

% Noise params down
window_duration_down = duration - start_down;
new_noise_down = noise_down(1:start_down/dt);

%% New simulations - UP
r_up_50 = zeros(1, duration);
r_up_neg50 = zeros(1, duration);

post_up = zeros(1, 21);
seed = 150; %20, 150
for freq = 1:21
    fprintf(' Current freq: %i\n', freq)
        
    if freq ~= 21
        stimuli(start_up/dt:start_up/dt + duration_s/dt - 1) = spike_train(freq,:);
        
        [window,window_t] = OUNoiseWindow(theta,0.5,window_duration_up,dt,save_dt,numsignals,noise_up(start_up/dt+1),seed);
        complete_window = [new_noise_up; window];


        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_up,x0,k,r0,stimuli), [0:dt:duration], ii_cc_up);
        r_s = ra_s(:,1)';
        a_s = ra_s(:,2)';
        [pre,post_stimuli] = durationMeasurement(r_s',start_up/dt, 'up');
        post_up(freq) = post_stimuli*dt;
        
        if freq == 5
            r_up_neg50 = r_s;
            stim_up_neg50 = stimuli;
        elseif freq == 15
            r_up_50 = r_s;
            stim_up_50 = stimuli;
        end
        
    end    
    if freq == 21
        stimuli(start_up/dt:start_up/dt + duration_s/dt - 1) = 0;
        
        [window,window_t] = OUNoiseWindow(theta,0.5,window_duration_up,dt,save_dt,numsignals,noise_up(start_up/dt+1),seed);
        complete_window = [new_noise_up; window];


        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_up,x0,k,r0,stimuli), [0:dt:duration], ii_cc_up);
        r_up_0 = ra_s(:,1)';
        a_s = ra_s(:,2)';
        [pre,post_nostimuli] = durationMeasurement(r_up_0',start_up/dt, 'up');
        post_up(freq) = post_nostimuli*dt;
        
        
    end    
end

% Rearrange
tmp = post_up;
post_up(11) = post_up(21);
post_up(12:end) = tmp(11:20);

% Fitting
f2 = -100:10:100;
lb = [0 -inf -100 0];
ub = [Inf Inf 100 Inf];
error_fun = @(p_new)(p_new(4) + (p_new(1)-p_new(4))./(1+exp(-p_new(2)*(f2-p_new(3)))))-post_up;
p0 = zeros(1,4);
p0(1) = max([post_up(1) post_up(end)]);
p0(4) = min([post_up(1) post_up(end)]);
slope = diff(post_up)/5;
[~, loc] = max(abs(slope));
p0(2) = slope(loc);
p0(3) = (f2(loc)+f2(loc+1))/2;

fitted_fun = @(p, x)(p(4) + (p(1)-p(4))./(1+exp(-p(2)*(x-p(3)))));

p_fitted = lsqnonlin(error_fun,p0,lb,ub);
fit_up = fitted_fun(p_fitted, [-100:0.1:100]);
sens_up = (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4));

figure; 
plot(f2, post_up, 'o'); hold on; plot([-100:0.1:100], fit_up);

%% New simulations - down
r_down_50 = zeros(1, duration);
r_down_neg50 = zeros(1, duration);
post_down = zeros(1, 21);
seed = 1; %1
for freq = 1:21
    fprintf(' Current freq: %i\n', freq)
        
    if freq ~= 21
        stimuli(start_down/dt:start_down/dt + duration_s/dt - 1) = spike_train(freq,:);
        
        [window ,window_t] = OUNoiseWindow(theta,0.5,window_duration_down,dt,save_dt,numsignals,noise_down(start_down/dt+1),seed);
        complete_window = [new_noise_down; window];


        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_down,x0,k,r0,stimuli), [0:dt:duration], ii_cc_down);
        r_s = ra_s(:,1)';
        a_s = ra_s(:,2)';
        [pre,post_stimuli] = durationMeasurement(r_s',start_down/dt, 'down');
        post_down(freq) = post_stimuli*dt;
        
        if freq == 5
            r_down_neg50 = r_s;
            stim_down_neg50 = stimuli;
        elseif freq == 15
            r_down_50 = r_s;
            stim_down_50 = stimuli;
        end
        
        
    end    
    if freq == 21
        stimuli(start_down/dt:start_down/dt + duration_s/dt - 1) = 0;
        
        [window,window_t] = OUNoiseWindow(theta,0.5,window_duration_down,dt,save_dt,numsignals,noise_down(start_down/dt+1),seed);
        complete_window = [new_noise_down; window];


        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_down,x0,k,r0,stimuli), [0:dt:duration], ii_cc_down);
        r_down_0 = ra_s(:,1)';
        a_s = ra_s(:,2)';
        [pre,post_nostimuli] = durationMeasurement(r_down_0',start_down/dt, 'down');
        post_down(freq) = post_nostimuli*dt;
        
    end    
end

% Rearrange
tmp = post_down;
post_down(11) = post_down(21);
post_down(12:end) = tmp(11:20);

% Fitting
f2 = -100:10:100;
lb = [0 -inf -100 0];
ub = [Inf Inf 100 Inf];
error_fun = @(p_new)(p_new(4) + (p_new(1)-p_new(4))./(1+exp(-p_new(2)*(f2-p_new(3)))))-post_down;
p0 = zeros(1,4);
p0(1) = max([post_down(1) post_down(end)]);
p0(4) = min([post_down(1) post_down(end)]);
slope = diff(post_down)/5;
[~, loc] = max(abs(slope));
p0(2) = slope(loc);
p0(3) = (f2(loc)+f2(loc+1))/2;

fitted_fun = @(p, x)(p(4) + (p(1)-p(4))./(1+exp(-p(2)*(x-p(3)))));

p_fitted = lsqnonlin(error_fun,p0,lb,ub);
fit_down = fitted_fun(p_fitted, [-100:0.1:100]);
sens_down = abs((p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4)));

figure; 
plot(f2, post_down, 'o'); hold on; plot([-100:0.1:100], fit_down);

%% Final plotting

f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17 20]*1.75);
fontSize = 16;
LineW = 2.5; 
stimuli = zeros(1,duration/dt+1); % Stimuli



% A1
ax1 = axes('Position',[0.11 0.89 0.375 0.08]);
ax1.PositionConstraint = 'innerposition';
plot(t,r_up_neg50,'LineWidth',LineW,'Color','k'); hold on; 
plot(t,stim_up_neg50,'LineWidth',LineW,'Color','#DAA520')
ylim([-1 1])
xlim([140 220])
xticklabels({''})
yticks([-1 0 1])
yticklabels({'-1 ','0 ','1 '})
ylabel('{\it r(t) / \color[rgb]{.85 .647 .125}s(t)}')
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

% A2
ax1 = axes('Position',[0.11 0.795 0.375 0.08]);
ax1.PositionConstraint = 'innerposition';
plot(t,r_up_0,'LineWidth',LineW,'Color','k'); hold on; 
plot(t,stimuli,'LineWidth',LineW,'Color','#DAA520')
ylim([-1 1])
xlim([140 220])
xticklabels({''})
yticks([-1 0 1])
yticklabels({'-1 ','0 ','1 '})
ylabel('{\it r(t) / \color[rgb]{.85 .647 .125}s(t)}')
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

% A3
ax1 = axes('Position',[0.11 0.7 0.375 0.08]);
ax1.PositionConstraint = 'innerposition';
plot(t,r_up_50,'LineWidth',LineW,'Color','k'); hold on; 
plot(t,stim_up_50,'LineWidth',LineW,'Color','#DAA520')
ylim([-1 1])
xlabel('Time (AU)')
xlim([140 220])
yticks([-1 0 1])
yticklabels({'-1 ','0 ','1 '})
ylabel('{\it r(t) / \color[rgb]{.85 .647 .125}s(t)}')
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

% A4
ax3 = axes('Position',[0.575 0.7 0.35 0.27]);
ax3.PositionConstraint = 'innerposition';
plot([-100:10:100],post_up,'.','MarkerSize',45,'Color','#027EDC','LineWidth',1); hold on; 
plot([-100:0.1:100],fit_up, 'Color', 'k', 'LineWidth', 2)
ylim([0 25])
ylabel('Duration (AU)')
xlabel('Frequency (Hz)')
title(sprintf('Sensibility: %.3f', sens_up))
set(ax3,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

%B1
ax1 = axes('Position',[0.11 0.515 0.375 0.08]);
ax1.PositionConstraint = 'innerposition';
plot(t,r_down_neg50,'LineWidth',LineW,'Color','k'); hold on; 
plot(t,stim_down_neg50,'LineWidth',LineW,'Color','#DAA520')
ylim([-1 1])
xlim([60 140])
xticklabels({''})
yticks([-1 0 1])
yticklabels({'-1 ','0 ','1 '})
ylabel('{\it r(t) / \color[rgb]{.85 .647 .125}s(t)}')
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

% B2
ax1 = axes('Position',[0.11 0.42 0.375 0.08]);
ax1.PositionConstraint = 'innerposition';
plot(t,r_down_0,'LineWidth',LineW,'Color','k'); hold on; 
plot(t,stimuli,'LineWidth',LineW,'Color','#DAA520')
ylim([-1 1])
xlim([60 140])
xticklabels({''})
yticks([-1 0 1])
yticklabels({'-1 ','0 ','1 '})
ylabel('{\it r(t) / \color[rgb]{.85 .647 .125}s(t)}')
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

% B3
ax1 = axes('Position',[0.11 0.325 0.375 0.08]);
ax1.PositionConstraint = 'innerposition';
plot(t,r_down_50,'LineWidth',LineW,'Color','k'); hold on; 
plot(t,stim_down_50,'LineWidth',LineW,'Color','#DAA520')
ylim([-1 1])
xlabel('Time (AU)')
xlim([60 140])
yticks([-1 0 1])
yticklabels({'-1 ','0 ','1 '})
ylabel('{\it r(t) / \color[rgb]{.85 .647 .125}s(t)}')
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')


% B4
ax4 = axes('Position',[0.575 0.325 0.35 0.27]);
ax4.PositionConstraint = 'innerposition';
plot([-100:10:100],post_down,'.','MarkerSize',45,'Color','#027EDC','LineWidth',1); hold on; 
plot([-100:0.1:100],fit_down, 'Color', 'k', 'LineWidth', 2)
ylim([0 25])
%yticks([0 5 10 15 20])
%yticklabels({'0','5','10','15','20'})
ylabel('Duration (AU)')
xlabel('Frequency (Hz)')
title(sprintf('Sensibility: %.3f', sens_down))
set(ax4,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

%% Plotting



%% Save figure ––––––– Uncomment and edit to save to personalised location

% cd '/Volumes/GoogleDrive-101271366273470520077/My Drive/PaperBelen/Figures/Temp figures'
% set(f,'Renderer','Painter')
% exportgraphics(gcf,'fig6.pdf','Resolution',300,'BackgroundColor','none')
% 



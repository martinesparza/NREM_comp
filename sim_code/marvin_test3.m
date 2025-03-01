%% Definitions.

%clear all
clc

% r = mean firing rate
% a = adaptation

% Arbitrary units for the time constants in order to avoid dimensions
tau_r = 1;
tau_a = 25; 

% Fixed values for sigmoidal function
x0 = 5; 
r0 = 0.5;
k = 15;

% Noise function is employed with Ornstein - Uhlenbeck method

duration = 600; %s
save_dt = 0.001; % Precision on saving values in noise function and in solving the ode
dt = 0.001; % must be the same as save_dt
theta =0.05; %About comparable...
sigma = 0.25;
numsignals = 1;

% load('noise_up.mat'); load('noise_t_up.mat'); load('ii_cc_up'); % UP
% load('noise_down.mat'); load('noise_t_down.mat'); load('ii_cc_down'); % DOWN

% load('spike_trains3.mat');
% load('spike_trains_10sec.mat');
% load('spike_trains_50sec.mat');



%% Bistable

% Weight values can change for different regimes
% Values for oscillatory regime
b = 1;
w = 6.3;
I = 2.35;

%% Oscillatory

% b = 1;
% w = 6;
% I = 2.5;

%% Inicializar el cluster

myCluster = parcluster('local');
myCluster.NumWorkers = 32;
parpool(32)

%% Simulacion

start = 100; %UP
% start = 50; %DOWN
noise_memory = zeros(600001,1);
new_noise = noise_up(1:start/dt); window_duration = duration - start;

duration_s = 10; %s
stimuli = zeros(1,600001);


post_durations = zeros(300,21);
avg_down = zeros(300,1);
avg_stimuli = zeros(300,1);


tic
for freq = 1:21
    
    post = zeros(200,1);
    
    if freq ~= 21
        stimuli(start/dt:start/dt + duration_s/dt - 1) = spike_trains_10sec(freq,:);
        
        parfor i = 1:300

        [window,window_t] = OUNoiseWindow(theta,0.5,window_duration,dt,save_dt,numsignals,noise_up(start/dt+1),29);
        complete_window = [new_noise; window];


        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_up,x0,k,r0,stimuli), [0:dt:duration], ii_cc);
        r_s = ra_s(:,1)';
        a_s = ra_s(:,2)';
        [pre,post_stimuli] = durationMeasurement(r_s',start/dt);
        post(i,1) = post_stimuli*dt;
        
        end
    end
    
    if freq == 21
        stimuli(start/dt:start/dt + duration_s/dt - 1) = 0;
        
        parfor i = 1:300

        [window,window_t] = OUNoiseWindow(theta,0.5,window_duration,dt,save_dt,numsignals,noise_down(start/dt+1),168);
        complete_window = [new_noise; window];


        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_down,x0,k,r0,stimuli), [0:dt:duration], ii_cc);
        r_s = ra_s(:,1)';
        a_s = ra_s(:,2)';
        [pre,post_nostimuli] = durationMeasurement(r_s',start/dt);
        post(i,1) = post_nostimuli*dt;
        
        avg_stimuli(i,1) = mean(complete_window(start/dt:start/dt+duration_s/dt,1));
        avg_down(i,1) = mean(complete_window(start/dt:start/dt+post_nostimuli,1));
        
        end
        
    end
    
     

    post_durations(:,freq) = post;
    
end

post_osc_down = zeros(300,21);
post_osc_down(:,1:10) = post_durations(:,1:10);
post_osc_down(:,11) = post_durations(:,21);
post_osc_down(:,12:21) = post_durations(:,11:20);

toc
save('wkspace_10_o_d.mat')

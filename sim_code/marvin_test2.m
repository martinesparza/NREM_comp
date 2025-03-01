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

load('noise_up.mat');
load('noise_t_up.mat');
load('ii_cc_up');
load('spike_trains.mat');


%% Bistable

% Weight values can change for different regimes
% Values for oscillatory regime
% b = 1;
% w = 6.3;
% I = 2.35;

%% Oscillatory

b = 1;
w = 6;
I = 2.5;

%% Inicializar el cluster

myCluster = parcluster('local');
myCluster.NumWorkers = 32;
parpool(32)

%% Simulacion

start = 100;
noise_memory = zeros(600001,1);
new_noise = noise_up(1:start/dt); window_duration = duration - start;

duration_s = 5; %s
stimuli = zeros(1,600001);


post_durations = zeros(200,42);


tic
for freq = 1:20
    stimuli(start/dt:start/dt + duration_s/dt - 1) = spike_trains(freq,:);
    post = zeros(200,1);
    
    parfor i = 1:200

        [window,window_t] = OUNoiseWindow(theta,0.5,window_duration,dt,save_dt,numsignals,noise_up(start/dt+1),i);
        complete_window = [new_noise; window];


        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_up,x0,k,r0,stimuli), [0:dt:duration], ii_cc);
        r_s = ra_s(:,1)';
        [pre,post_stimuli] = durationMeasurement(r_s',start/dt);

        post(i,1) = post_stimuli*dt;
        

    end 

    post_durations(:,freq) = post;
    
end

toc
save('workspace2.mat')
fprintf('Prueba 2 completada.\n')




%%
% clear all
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

duration = 5000; %s
save_dt = 0.01; % Precision on saving values in noise function and in solving the ode
dt = 0.01; % must be the same as save_dt
theta =0.05; %About comparable...
sigma = 0.0182;
numsignals = 1;

% [noise_down, noise_t_down] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
% load('noise_up.mat');
% load('noise_t_up.mat');
% load('ii_cc_up.mat');

load('noise_down.mat');
load('noise_t_down.mat');
load('ii_cc_down.mat');

load('spike_train.mat');

%% Inicializar el cluster
delete(gcp('nocreate'))
myCluster = parcluster('local');
myCluster.NumWorkers = 32;
parpool(32)

%% Simulacion

% start = 1020; % UP 
start = 600; %DOWN 
ii_cc = ii_cc_down;
new_noise = noise_down(1:start/dt); window_duration = duration - start;

duration_s = 10; %s
stimuli = zeros(1,500001);

post_down2_bis = zeros(300,21);
avg_down2_bis = zeros(300,1);
post_down2_osc = zeros(300,21);
avg_down2_osc = zeros(300,1);

%% Bis
%% Bistable

b = 1;
w = 6.3;
I = 2.35;

post_durations = zeros(300,21);

for freq = 1:21
    if freq ~= 21
        
        stimuli(start/dt:start/dt + duration_s/dt - 1) = spike_train(freq,:);
        post = zeros(300,1);
        avg = zeros(300,1);
        
        parfor i = 1:300

        [window,window_t] = OUNoiseWindow(theta,sigma,window_duration,dt,save_dt,numsignals,noise_down(start/dt+1),i);
        complete_window = [new_noise; window];


        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_down,x0,k,r0,stimuli), [0:dt:duration], ii_cc);
        r_s = ra_s(:,1)';
        a_s = ra_s(:,2)';
        [pre,post_stimuli] = durationMeasurement(r_s',start/dt);
        post(i,1) = post_stimuli*dt;
        
        end
    end
    
    if freq == 21
        stimuli(start/dt:start/dt + duration_s/dt - 1) = 0;
        post = zeros(300,1);
        avg = zeros(300,1);
        
        parfor i = 1:300

        [window,window_t] = OUNoiseWindow(theta,sigma,window_duration,dt,save_dt,numsignals,noise_down(start/dt+1),i);
        complete_window = [new_noise; window];


        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_down,x0,k,r0,stimuli), [0:dt:duration], ii_cc);
        r_s = ra_s(:,1)';
        a_s = ra_s(:,2)';
        [pre,post_nostimuli] = durationMeasurement(r_s',start/dt);
        post(i,1) = post_nostimuli*dt;
        
%         avg(i,1) = mean(complete_window(start/dt:start/dt+duration_s/dt,1));
        avg(i,1) = mean(complete_window(start/dt:start/dt+post_nostimuli,1));
        
        end
        
    end
    
     

    post_durations(:,freq) = post;
    
end

post_down2_bis(:,1:10) = post_durations(:,1:10);
post_down2_bis(:,11) = post_durations(:,21);
post_down2_bis(:,12:21) = post_durations(:,11:20);

avg_down2_bis(:,1) = avg;

%% Osc
%% Oscillatory

b = 1;
w = 6;
I = 2.5;

post_durations = zeros(300,21);

for freq = 1:21
    
    if freq ~= 21
        
        stimuli(start/dt:start/dt + duration_s/dt - 1) = spike_train(freq,:);
        post = zeros(300,1);
        avg = zeros(300,1);
        
        parfor i = 1:300

        [window,window_t] = OUNoiseWindow(theta,sigma,window_duration,dt,save_dt,numsignals,noise_down(start/dt+1),i);
        complete_window = [new_noise; window];


        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_down,x0,k,r0,stimuli), [0:dt:duration], ii_cc);
        r_s = ra_s(:,1)';
        a_s = ra_s(:,2)';
        [pre,post_stimuli] = durationMeasurement(r_s',start/dt);
        post(i,1) = post_stimuli*dt;
        
        end
    end
    
    if freq == 21
        stimuli(start/dt:start/dt + duration_s/dt - 1) = 0;
        post = zeros(300,1);
        avg = zeros(300,1);
        
        parfor i = 1:300

        [window,window_t] = OUNoiseWindow(theta,sigma,window_duration,dt,save_dt,numsignals,noise_down(start/dt+1),i);
        complete_window = [new_noise; window];


        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_down,x0,k,r0,stimuli), [0:dt:duration], ii_cc);
        r_s = ra_s(:,1)';
        a_s = ra_s(:,2)';
        [pre,post_nostimuli] = durationMeasurement(r_s',start/dt);
        post(i,1) = post_nostimuli*dt;
        
%         avg(i,1) = mean(complete_window(start/dt:start/dt+duration_s/dt,1));
        avg(i,1) = mean(complete_window(start/dt:start/dt+post_nostimuli,1));
        
        end
        
    end
    
     

    post_durations(:,freq) = post;
    
end


post_down2_osc(:,1:10) = post_durations(:,1:10);
post_down2_osc(:,11) = post_durations(:,21);
post_down2_osc(:,12:21) = post_durations(:,11:20);

avg_down2_osc(:,1) = avg;

%% Save variables

save('avg_down2_bis.mat','avg_down2_bis')
save('post_down2_bis.mat','post_down2_bis')
save('avg_down2_osc.mat','avg_down2_osc')
save('post_down2_osc.mat','post_down2_osc')

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

duration = 5000; %s
save_dt = 0.01; % Precision on saving values in noise function and in solving the ode
dt = 0.01; % must be the same as save_dt
theta =0.05; %About comparable...
% sigma = 0.075;
% sigma = 0.25;
sigma = 0.0096;
numsignals = 1;

% [noise_up, noise_t_up] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
load('noise_up.mat');
load('noise_t_up.mat');
load('ii_cc_up.mat');

% load('noise_down.mat');
% load('noise_t_down.mat');
% load('ii_cc_down.mat');

%% Inicializar el cluster
delete(gcp('nocreate'))
myCluster = parcluster('local');
myCluster.NumWorkers = 32;
parpool(32)

%% Simulacion

start = 1020; %UP 
% start = 200; %DOWN 0.07
% ii_cc = ii_cc_down;
new_noise = noise_up(1:start/dt); window_duration = duration - start;

post_up_bis = zeros(300,1);
avg_up_bis = zeros(300,1);
post_up_osc = zeros(300,1);
avg_up_osc = zeros(300,1);


%% Bis
%% Bistable

% Weight values can change for different regimes
% Values for oscillatory regime
b = 1;
w = 6.3;
I = 2.35;

post = zeros(300,1);
avg = zeros(300,1);
    
parfor i = 1:300

[window,window_t] = OUNoiseWindow(theta,sigma,window_duration,dt,save_dt,numsignals,noise_up(start/dt+1),i);
complete_window = [new_noise; window];

noise_t_up = [0:dt:duration];
[t,ra_s] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_down,x0,k,r0), [0:dt:duration], ii_cc);
r_s = ra_s(:,1)';
a_s = ra_s(:,2)';
[pre,post_nostimuli] = durationMeasurement(r_s',start/dt);
post(i,1) = post_nostimuli*dt;
avg(i,1) = mean(complete_window(start/dt:start/dt+post_nostimuli,1));

end

post_up_bis(:,1) = post;
avg_up_bis(:,1) = avg;

%% Osc

%% Oscillatory

b = 1;
w = 6;
I = 2.5;


post = zeros(300,1);
avg = zeros(300,1);
    
clear window window_t complete_window
parfor i = 1:300

[window,window_t] = OUNoiseWindow(theta,sigma,window_duration,dt,save_dt,numsignals,noise_up(start/dt+1),i);
complete_window = [new_noise; window];

noise_t_up = [0:dt:duration];
[t,ra_s] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,complete_window,noise_t_up,x0,k,r0), [0:dt:duration], ii_cc);
r_s = ra_s(:,1)';
a_s = ra_s(:,2)';
[pre,post_nostimuli] = durationMeasurement(r_s',start/dt);
post(i,1) = post_nostimuli*dt;
avg(i,1) = mean(complete_window(start/dt:start/dt+post_nostimuli,1));

end


post_up_osc(:,1) = post;
avg_up_osc(:,1) = avg;    


save('avg_up_bis.mat','avg_up_bis')
save('post_up_bis.mat','post_up_bis')
save('avg_up_osc.mat','avg_up_osc')
save('post_up_osc.mat','post_up_osc')

%% PLOT__ COMMENT FOR MARVIN

% figure; 
% plot(avg_up_bis,post_up_bis,'bo','LineWidth',2.0,'MarkerSize',12);
% xlabel('Noise'); ylabel('Durations');
% title('UP. BISTABLE');
% xlim([-0.04 0.04])
% ylim([1 10000])
% set(gca,'yscale','log')
% set(gca,'FontSize',24,'box','off','LineWidth',2)
% set(gcf,'Position',[0 0 800 600])
% % 
% figure; 
% plot(avg_up_osc,post_up_osc,'bo','LineWidth',2.0,'MarkerSize',12);
% xlabel('Noise'); ylabel('Durations');
% title('UP. OSCILLATORY');
% xlim([-0.04 0.04])
% ylim([1 1000])
% set(gca,'yscale','log')
% set(gca,'FontSize',24,'box','off','LineWidth',2)
% set(gcf,'Position',[0 0 800 600])
% 
% figure; 
% plot(avg_down_bis,post_down_bis,'ro','LineWidth',2.0,'MarkerSize',12);
% xlabel('Noise'); ylabel('Durations');
% title('DOWN. BISTABLE');
% xlim([-0.6 0.8])
% ylim([1 1000])
% set(gca,'yscale','log')
% set(gca,'FontSize',24,'box','off','LineWidth',2)
% set(gcf,'Position',[0 0 800 600])

% figure; 
% plot(avg_down_osc,post_down_osc,'ro','LineWidth',2.0,'MarkerSize',12); hold on;
% % plot(avg_25_osc,post_25_osc,'ko','LineWidth',2.0,'MarkerSize',12,'MarkerFaceColor','k');
% xlabel('Noise'); ylabel('Durations');
% title('DOWN. OSCILLATORY');
% % xlim([-0.6 0.8])
% ylim([1 1000])
% % legend('sigma = 0.07','sigma = 0.25')
% set(gca,'yscale','log')
% set(gca,'FontSize',24,'box','off','LineWidth',2)
% set(gcf,'Position',[0 0 800 600])


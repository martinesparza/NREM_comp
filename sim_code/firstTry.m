%% Definitions
clear all
clc
close all
% r = mean firing rate
% a = adaptation

% Arbitrary units for the time constants in order to avoid dimensions
tau_r = 1;
tau_a = 25; 

% Fixed values for sigmoidal function
x0 = 5; 
r0 = 0.5;
k = 15;

% Weight values can change for different regimes
% Values for oscillatory regime
b = 1;
w = 6.3;
I = 2.35;

% Noise function is employed with Ornstein - Uhlenbeck method

duration = 1000; %s
save_dt = 0.01; % Precision on saving values in noise function and in solving the ode
dt = 0.01; % must be the same as save_dt
theta =0.05; %About comparable...
sigma = 0.25;
numsignals = 1;

[noise, noise_t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals,1);

% New input. Stimuli

amplitude = 1;
duration_s = 1000; %10s 
% start = randi(100000 - duration);
% interval = start:(start + duration);

stimuli = zeros(1,100001);
stimuli(60000:70000) = amplitude; %Amplitude of stimul

% Initial conditions



%% System
tic
[t,ra] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0), [0:dt:duration], [0.5; 0]);
r = ra(:,1);
a = ra(:,2);

figure
subplot(2,1,1)
plot(t,r,'-o')
title('Time course. Firing rate.')
xlabel('Time (AU)')
ylabel('Rate')

subplot(2,1,2)
plot(t,a,'-o')
title('Time course. Adaptation.')
xlabel('Time (AU)')
ylabel('Adaptation')

toc

%% System with stimuli

tic
[t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0,stimuli), [0:dt:duration], [0.5; 0]);
r_s = ra_s(:,1);
a_s = ra_s(:,2);

figure
subplot(2,1,1)
plot(t,r_s,'-o')
title('Time course. Firing rate.')
xlabel('Time (AU)')
ylabel('Rate')

subplot(2,1,2)
plot(t,a_s,'-o')
title('Time course. Adaptation.')
xlabel('Time (AU)')
ylabel('Adaptation')

toc
%% Nullclines

sample_r = 0:0.001:1; % Sample r vector to draw nullclines. 


nul_a = 1./(1 + exp(-k.*(sample_r-r0)));
nul_r = (1/b) * (w.*sample_r + I - x0 + log(((1./sample_r) - 1)));



figure
plot(sample_r,nul_a,'b',sample_r,nul_r,'k','LineWidth',1.5); hold on;
% plot(sample_r(locs(1,1)),nul_a(locs(1,1)),'wo','LineWidth',1.5,'MarkerFaceColor','w','MarkerEdgeColor','k')
axis([0 1 0 1])
legend('da/dt = 0','dr/dt = 0')
xlabel('Pop. rate, r')
ylabel('Adaptation, a')
title('Phase plane')


%% Thresholding

tic
[r_thresh,r_cross,r_bihist,r_diptest] = BimodalThresh(r,'Schimdt'); % Calculate treshold and up-down transitions. 
toc

tic
[a_thresh,a_cross,a_bihist,a_diptest] = BimodalThresh(a,'Schmidt'); % Calculate treshold and up-down transitions. 
toc
%% New plotting

% Plotting over the temporal solutions, the up and down intervals

dim_up_r = length(r_cross.upints);
dim_down_r = length(r_cross.downints);

dim_up_a = length(a_cross.upints);
dim_down_a = length(a_cross.downints);

figure
subplot(2,1,1)
plot(t,r,'-o'); hold on;
for i = 1:dim_up_r
plot(t(r_cross.upints(i,1):r_cross.upints(i,2)),r(r_cross.upints(i,1):r_cross.upints(i,2)),'-bo'); hold on
end
for i = 1:dim_down_r
plot(t(r_cross.downints(i,1):r_cross.downints(i,2)),r(r_cross.downints(i,1):r_cross.downints(i,2)),'-ro'); hold on
end

subplot(2,1,2)
plot(t,a,'-o'); hold on;
for i = 1:dim_up_a
plot(t(a_cross.upints(i,1):a_cross.upints(i,2)),a(a_cross.upints(i,1):a_cross.upints(i,2)),'-bo'); hold on
end
for i = 1:dim_down_a
plot(t(a_cross.downints(i,1):a_cross.downints(i,2)),a(a_cross.downints(i,1):a_cross.downints(i,2)),'-ro'); hold on
end


%% 100 simulations

% Parameters of stimuli

%amplitude = 1;
duration_s = 10000; %100s 
stimuli = zeros(1,100001);
start = 60000;
finish = start+duration_s;
stimuli(start:finish) = amplitude; 
margin = duration_s;

% memory_up = zeros(10,(duration_s + (2*margin) + 1)); % memory matrix to save data
% memory_down = zeros(10,20001);
n_up = 0; % counter for every row.
n_down = 0;

%% Initialize cell memory

start = 60000;
stimuli = zeros(1,100001);

max_sim = 4; % total number of simulation (100) 
amplitude_vec = 0.5:0.1:0.5;
duration_vec = 10:10:10; % This means that duration will be simulated from 10s to 100s in intervals of 10s
memory_cell_up = initialize('Stimulation time (s) = ',duration_vec,amplitude_vec,max_sim,dt);

 % counter for every row.
n_down = 0;

% While loop for temporal average. It can be for "r" or for "a".

% UP CASE
tic

for j = 1:length(duration_vec)
duration_s = duration_vec(j)/dt;
margin = 100;
margin_i = margin/dt;


    for i = 1:length(amplitude_vec)
        stimuli(start:start + duration_s) = amplitude_vec(i);

        parfor parall = 1:max_sim
            n_up = 0;
            while n_up < 1
                [noise, noise_t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
                ic = rand(2,1);

                [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0,stimuli), [0:dt:duration], ic);
                r_s = ra_s(:,1);
                %a_s = ra_s(:,2);

                [rs_thresh,rs_cross,rs_bihist,rs_diptest] = BimodalThresh(r_s,'Schimdt');
                %[a_thresh,a_cross,a_bihist,a_diptest] = BimodalThresh(a_s,'Schmidt');

                dim_up_r = length(rs_cross.upints);
                N = histcounts(start,reshape(rs_cross.upints',[],1)); % Vamos a comprobar si el inicio esta en un intervalo up

                if mod(find(N == 1),2) ~= 0
                    n_up = 1
                    memory_cell_up{j,2}(parall,:,i) = r_s((start-margin_i):(start + duration_s +margin_i));
                end
            end
        end
    end
end
toc

% aux_t = -margin:(duration_s+margin); %Auxiliary temporal vector 
% figure;
% plot(aux_t,memory_up(1:n_up,:)');
%% DOWN CASE

while n_down < 3
    [noise, noise_t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
    ic = rand(2,1);

    [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0,stimuli), [0:dt:duration], ic);
    r_s = ra_s(:,1);
    %a_s = ra_s(:,2);

    [rs_thresh,rs_cross,rs_bihist,rs_diptest] = BimodalThresh(r_s,'Schimdt');
    %[a_thresh,a_cross,a_bihist,a_diptest] = BimodalThresh(a_s,'Schmidt');

    dim_down_r = length(rs_cross.downints); d = false; 

    for i=1:dim_down_r
        if (rs_cross.downints(i,1) <= start) && (rs_cross.downints(i,2) >= start)
            d = true; % Our impulse starts in a down interval
            n_down = n_down + 1;
        end
    end

    if d
        memory_down(n_down,:) = r_s(55000:75000);
    end
end


%% PDF STUDY



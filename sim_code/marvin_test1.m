%% Definitions

clear all
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

load('noise_marvin.mat');
load('noise_t_marvin.mat');

% [noise, noise_t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);

%%
ic = [0.5;0.5];
[t,ra_s] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0), [0:dt:duration], ic);
r_s = ra_s(:,1);

start = 200/dt;
[ pre,post ] = durationMeasurement( r_s,start );
final = pre + post;
save('workspace.mat')
fprintf('Prueba 1 completada.\n')




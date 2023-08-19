function [ dradt,r_inf  ] = rateAdapatationStimuli( t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0,stimuli)
% Funcion de sistema de ecuaciones diferenciales para resolver mediante
% ode45. Ademas en este caso anyadiremos un estimulo que iremos variando. 
%   Estaremos integrando la mean firing rate de la poblacion y la
%   adaptacion del sistema

noise_inter = interp1(noise_t,noise,t); % Interpol of noise 
stimuli_inter = interp1(noise_t,stimuli,t); % Interpol of stimuli

r_inf = (1/(1 + exp(-(w*ra(1) - b*ra(2) + (I + stimuli_inter) + noise_inter - x0))));

dradt = zeros(2,1);
dradt(1) = (1/tau_r) * ( -ra(1) + r_inf);
dradt(2) = (1/tau_a) * ( -ra(2) + (1/(1 + exp(-k*(ra(1) - r0)))));



end


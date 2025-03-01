function status = myOutputFcnStimuli(t,ra,flag,w,b,I,noise,noise_t,x0,stimuli)
% OutputFcn sample
noise_inter = interp1(noise_t,noise,t);
stimuli_inter = interp1(noise_t,stimuli,t); % Interpol of stimuli

persistent r_inf

switch flag
    case 'init'
        r_inf = (1./(1 + exp(-(w*ra(1) - b*ra(2) + (I + stimuli_inter) + noise_inter - x0)))); %mass increases linearly with time
    case ''
        r_inf = [r_inf, (1./(1 + exp(-(w*ra(1) - b*ra(2) + (I + stimuli_inter) + noise_inter - x0))))];
    case 'done'
        assignin('base','r_inf',r_inf);
end

status = 0;


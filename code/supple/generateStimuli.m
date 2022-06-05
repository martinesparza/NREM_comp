function [spike_trains,spikeMat,tVec] = generateStimuli(duration_s, dt,trials, freqs)
%generateStimuli Generate the stimuli for the model using a spike train and
%a convolution with the alpha function
%   Detailed explanation goes here

positive_currents = zeros(length(freqs),duration_s/dt);

for i = 1:length(freqs)
    [spikeMat, tVec] = poissonSpikeGen(freqs(i), duration_s, trials, dt);
    spikeMat2 = double(spikeMat);
    alfa = 1/0.1;
    alfa_function = alfa^2 .* tVec(1:250) .* exp(-alfa.*tVec(1:250));
    positive_current = conv(spikeMat2,alfa_function);
    positive_current = positive_current(1:duration_s/dt)./max(positive_current);
    positive_currents(i,:) = positive_current;
end

negative_currents = -positive_currents;
spike_trains = [flip(negative_currents,1); positive_currents];


end


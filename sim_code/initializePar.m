function [out] = initializePar(duration,amplitude_vec,max_sim,dt)
% Function intended to initialize a series of matrices
%   Detailed explanation goes here
z_axis = length(amplitude_vec);
iterations = duration/dt;
margin = 10000;

out = zeros(max_sim, iterations + (2*margin) + 1,z_axis);

end


function dradt = model(t,ra,mu,noise)
% Funcion de sistema de ecuaciones diferenciales para resolver mediante
% ode45
%   Detailed explanation goes here

sample_t = round(t);
dradt = zeros(2,1);
dradt(1) = ra(2) + noise(sample_t,:);
dradt(2) = mu*(1-ra(1)^2)*ra(2)-ra(1);
dradt;
end


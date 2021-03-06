function [ X,T,aux,randnums ] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals,seed)
%[ X,T ] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals)) simulates
%Ornstein-Uhlenbeck noise with standard deviation sigma and characteristic
%frequency theta.
%Initial conditions: random value from normal distribution around 0 with
%std=sigma
%Simulated for a given duration with timestep dt, but output is only saved 
%at time points save_dt. Can simulate multiple signals simultaneously.
%
%
%NOTE: variance is given by sigma^2/(2*theta)
%DLevenstein March2018
%TO DO: Make save_dt, numsignals optional
%% DEV
% duration = 20; %s
% save_dt = 1;
% dt = 0.01;
% theta =0.05; %About comparable...
% sigma = 0.25;
% numsignals = 1;



%%
rng(seed)
simtimevector = 0:dt:duration;
SimTimeLength = length(simtimevector);
randnums = randn(numsignals,SimTimeLength);
savecounter = 1;
aux = randn(numsignals,1);
X_t = sigma.*aux; %Start at random value - gaussian distribted around sigma
clear X
clear T
for tt=1:SimTimeLength
    dX = -theta.*X_t.*dt + sqrt(2.*theta).*sigma.*randnums(:,tt).*sqrt(dt);
    %dX = sqrt(dt)*sigma.*randnums(:,tt); % White noise. 
    X_t = X_t + dX;
    %X_t = dX;
    
    if mod(simtimevector(tt),save_dt)==0
        X(savecounter,:) = X_t; %Pre-allocate
        T(savecounter,:) = simtimevector(tt);
        savecounter = savecounter+1;
    end
end


% figure
% plot(T,X)
%add figure of mean/std/power spectrum. with SHOWFIG=TRUE

end


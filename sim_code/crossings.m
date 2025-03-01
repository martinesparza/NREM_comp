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


sigmas = linspace(0.05,0.25,10);

% Noise function is employed with Ornstein - Uhlenbeck method

duration = 50000; %s
save_dt = 0.01; % Precision on saving values in noise function and in solving the ode
dt = 0.01; % must be the same as save_dt
 %About comparable...
theta = 0.05;

numsignals = 1;

up_sigmas_bis = cell(1,2);
down_sigmas_bis = cell(1,2);

% for i = 1:10

% sigma = sigmas(i);
sigma = 0.25;



% load('ii_cc_down');

%% Bistable

% Weight values can change for different regimes
% Values for oscillatory regime
b = 1;
w = 6.3;
I = 2.35;

%% Oscillatory

% b = 1;
% w = 6;
% I = 2.5;

%% Inicializar el cluster
delete(gcp('nocreate'))
myCluster = parcluster('local');
myCluster.NumWorkers = 32;
parpool(32)

%% Simulacion
tic

upcrossing_bis = 0;
downcrossing_bis = 0;
% parfor j = 1:50

[noise, noise_t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
rng(1)
[t,ra_s] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0), [0:dt:duration], rand(2,1));
r_s = ra_s(:,1)';
save('r.mat','r_s')

[rs_thresh,rs_cross,rs_bihist,rs_diptest] = BimodalThresh(r_s,'Schimdt');

upcrossing_bis = [upcrossing_bis; noise(rs_cross.upints(1,:))];
downcrossing_bis = [downcrossing_bis; noise(rs_cross.upints(2,:))];

% end

up_sigmas_bis{1,1} = upcrossing_bis;
down_sigmas_bis{1,1} = downcrossing_bis;
toc
% end


save('up_sigmas_bis.mat','up_sigmas_bis')
save('down_sigmas_bis.mat','down_sigmas_bis')
save('r.mat','r_s')

%% Plot. COMENTAR ESTA SECCION PARA MARVIN
sigmas = linspace(0.05,0.25,10);
% % sigma = 
% % Plot hist
% 
% for i = 1
%     figure;
%     histogram(up_sigmas_osc{1,i},'Normalization','probability'); hold on; histogram(down_sigmas_osc{1,i},'Normalization','probability');
%     xlabel('Time')
%     ylim([0 0.2])
%     xlim([-1 1])
%     legend('UP','DOWN')
%     xlabel('Noise')
%     str = sprintf('Oscillatory. Sigma: %.2f',sigmas(i));
%     title(str)
%     set(gca,'FontSize',24,'box','off','LineWidth',2)
%     set(gcf,'Position',[0 0 800 600])
%     saveas(gcf,strcat(str,'.png'))
% end

% Plot Mean & std
means_osc = zeros(2,10);
stds_osc = zeros(2,10);

means_bis = zeros(2,10);
stds_bis = zeros(2,10);

for i = 1:10
    means_osc(1,i) = mean(up_sigmas_osc{1,i});
    stds_osc(1,i) = std(up_sigmas_osc{1,i}); 
    means_osc(2,i) = mean(down_sigmas_osc{1,i});
    stds_osc(2,i) = std(down_sigmas_osc{1,i});
    
    means_bis(1,i) = mean(up_sigmas_bis{1,i});
    stds_bis(1,i) = std(up_sigmas_bis{1,i}); 
    means_bis(2,i) = mean(down_sigmas_bis{1,i});
    stds_bis(2,i) = std(down_sigmas_bis{1,i});
end

figure;
plot(sigmas,means_osc(1,:),'-b','LineWidth',2.0); hold on;
plot(sigmas,means_osc(2,:),'-r','LineWidth',2.0); hold on; 
plot(sigmas,means_bis(1,:),'--b','LineWidth',2.0); hold on;
plot(sigmas,means_bis(2,:),'--r','LineWidth',2.0); hold on;
xlabel('Sigma')
legend('Osc. UP','Osc. DOWN','Bis. UP','Bis. DOWN')
ylabel('Noise')
title('Mean vs. Sigma')
set(gca,'FontSize',24,'box','off','LineWidth',2)
set(gcf,'Position',[0 0 800 600])
% saveas(gcf,'Mean.png')

figure;
plot(sigmas,stds_osc(1,:),'-b','LineWidth',2.0); hold on;
plot(sigmas,stds_osc(2,:),'-r','LineWidth',2.0); hold on; 
plot(sigmas,stds_bis(1,:),'--b','LineWidth',2.0); hold on;
plot(sigmas,stds_bis(2,:),'--r','LineWidth',2.0); hold on;
xlabel('Sigma')
legend('Osc. UP','Osc. DOWN','Bis. UP','Bis. DOWN')
ylabel('Noise')
title('STD vs. Sigma')
set(gca,'FontSize',24,'box','off','LineWidth',2)
set(gcf,'Position',[0 0 800 600])
% saveas(gcf,'STD.png')

figure;
plot(sigmas,stds_osc(1,:)./means_osc(1,:),'-b','LineWidth',2.0); hold on;
plot(sigmas,stds_osc(2,:)./means_osc(2,:),'-r','LineWidth',2.0); hold on; 
plot(sigmas,stds_bis(1,:)./means_bis(1,:),'--b','LineWidth',2.0); hold on;
plot(sigmas,stds_bis(2,:)./means_bis(2,:),'--r','LineWidth',2.0); hold on;
xlabel('Sigma')
legend('Osc. UP','Osc. DOWN','Bis. UP','Bis. DOWN')
ylabel('Noise')
title('Variation Coefficient vs. Sigma')
set(gca,'FontSize',24,'box','off','LineWidth',2)
set(gcf,'Position',[0 0 800 600])
saveas(gcf,'CV.png')





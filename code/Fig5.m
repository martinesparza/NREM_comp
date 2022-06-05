%% Figures @NREM_comp
%
% Theoretical and Computational Neuroscience (TCN), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 5. Differences in noise amplitude reveal the influence of the noise 
% term on r(t) in the bistable and oscillatory regimes
%
% Created: Mon 1 Mar 2021, 09:28
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Mon 23 May 2022, 17:22
% Last edited by: Martin Esparza-Iaizzo

%% Load data

% Point to data location
addpath(genpath('/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/final_data_00182/UP'))
% Load data
tmp = load('avg_up2_osc.mat');
sigma_low.avg_up_osc = tmp.avg_up_osc;
tmp = load('avg_up2_bis.mat');
sigma_low.avg_up_bis = tmp.avg_up_bis;
tmp = load('post_up2_bis.mat');
sigma_low.post_up_bis = tmp.post_up_bis;
tmp = load('post_up2_osc.mat');
sigma_low.post_up_osc = tmp.post_up_osc;


addpath(genpath('/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/final_data_00182/DOWN'))
tmp = load('avg_down2_osc.mat');
sigma_low.avg_down_osc = tmp.avg_down_osc;
tmp = load('avg_down2_bis.mat');
sigma_low.avg_down_bis = tmp.avg_down_bis;
tmp = load('post_down2_bis.mat');
sigma_low.post_down_bis = tmp.post_down_bis;
tmp = load('post_down2_osc.mat');
sigma_low.post_down_osc = tmp.post_down_osc;

addpath(genpath('/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/sigma = 0.07/0.25/UP - 0.25'))
tmp = load('avg_up_osc.mat');
sigma_high.avg_up_osc = tmp.avg_up_osc;
tmp = load('avg_up_bis.mat');
sigma_high.avg_up_bis = tmp.avg_up_bis;
tmp = load('post_up_bis.mat');
sigma_high.post_up_bis = tmp.post_up_bis;
tmp = load('post_up_osc.mat');
sigma_high.post_up_osc = tmp.post_up_osc;

addpath(genpath('/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/sigma = 0.07/0.25/DOWN - 0.25'))
tmp = load('avg_down_osc.mat');
sigma_high.avg_down_osc = tmp.avg_down_osc;
tmp = load('avg_down_bis.mat');
sigma_high.avg_down_bis = tmp.avg_down_bis;
tmp = load('post_down_bis.mat');
sigma_high.post_down_bis = tmp.post_down_bis;
tmp = load('post_down_osc.mat');
sigma_high.post_down_osc = tmp.post_down_osc;

%% Plot

f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17 20]*1.75);
fontSize = 16;
MarkerSize = 10;

% A
ax1 = axes('Position',[0.15 0.67 0.36 0.25]);
ax1.PositionConstraint = 'innerposition';
plot(sigma_high.avg_down_bis,sigma_high.post_down_bis,'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on; 
plot(sigma_high.avg_down_osc,sigma_high.post_down_osc,'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8')
ylim([1 10000])
xlim([-0.5 0.5])
xlabel('<\xi>_t')
ylabel('Down duration (AU)')
yticklabels({'10^0','','10^2','','10^4'})
yticks([1 10 100 1000 10000])
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')

% B
ax2 = axes('Position',[0.575 0.67 0.36 0.25]);
ax2.PositionConstraint = 'innerposition';
plot(sigma_high.avg_up_bis,sigma_high.post_up_bis,'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on; 
plot(sigma_high.avg_up_osc,sigma_high.post_up_osc,'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8')
ylim([1 10000])
xlim([-0.5 0.5])
xlabel('<\xi>_t')
ylabel('Up Duration (AU)')
yticklabels({''})
yticks([1 10 100 1000 10000])
set(ax2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')

% C
ax3 = axes('Position',[0.15 0.3 0.36 0.25]);
ax3.PositionConstraint = 'innerposition';
plot(sigma_low.avg_down_bis,sigma_low.post_down_bis(:,11),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on; 
plot(sigma_low.avg_down_osc,sigma_low.post_down_osc(:,11),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8')
ylim([1 10000])
xlim([-0.05 0.05])
xlabel('<\xi>_t')
ylabel('Down duration (AU)')
yticklabels({'10^0','','10^2','','10^4'})
yticks([1 10 100 1000 10000])
set(ax3,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')


% D
ax4 = axes('Position',[0.575 0.3 0.36 0.25]);
ax4.PositionConstraint = 'innerposition';
plot(sigma_low.avg_up_bis,sigma_low.post_up_bis(:,11),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on; 
plot(sigma_low.avg_up_osc,sigma_low.post_up_osc(:,11),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8')
ylim([1 10000])
xlim([-0.05 0.05])
xlabel('<\xi>_t')
ylabel('Up Duration (AU)')
yticklabels({''})
yticks([1 10 100 1000 10000])
set(ax4,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')


%% Save figure ––––––– Uncomment and edit to save to personalised location

% cd '/Volumes/GoogleDrive-101271366273470520077/My Drive/PaperBelen/Figures/Temp figures'
% set(f,'Renderer','Painter')
% exportgraphics(gcf,'fig5.pdf','Resolution',300,'BackgroundColor','none')


%% Stats –– Pearson correlation

% Clear variables
clear rho_sigma1_DOWN pval_sigma1_DOWN rho_sigma1_UP pval_sigma1_UP
clear rho_sigma2_DOWN pval_sigma2_DOWN rho_sigma2_UP pval_sigma2_UP


% Sigma 1 – DOWN
[rho_sigma1_DOWN.bis,pval_sigma1_DOWN.bis] = corr(sigma_high.avg_down_bis,sigma_high.post_down_bis); % Bistable
[rho_sigma1_DOWN.osc,pval_sigma1_DOWN.osc] = corr(sigma_high.avg_down_osc,sigma_high.post_down_osc); % Oscillatory

% Sigma 1 - UP
[rho_sigma1_UP.bis,pval_sigma1_UP.bis] = corr(sigma_high.avg_up_bis,sigma_high.post_up_bis); % Bistable
[rho_sigma1_UP.osc,pval_sigma1_UP.osc] = corr(sigma_high.avg_up_osc,sigma_high.post_up_osc); % Oscillatory

% Sigma 2 – DOWN
[rho_sigma2_DOWN.bis,pval_sigma2_DOWN.bis] = corr(sigma_low.avg_down_bis,sigma_low.post_down_bis(:,11)); % Bistable
[rho_sigma2_DOWN.osc,pval_sigma2_DOWN.osc] = corr(sigma_low.avg_down_osc,sigma_low.post_down_osc(:,11)); % Oscillatory

% Sigma 2 - UP
[rho_sigma2_UP.bis,pval_sigma2_UP.bis] = corr(sigma_low.avg_up_bis,sigma_low.post_up_bis(:,11)); % Bistable
[rho_sigma2_UP.osc,pval_sigma2_UP.osc] = corr(sigma_low.avg_up_osc,sigma_low.post_up_osc(:,11)); % Oscillatory



%% Figures @NREM_comp
%
% Theoretical and Computational Neuroscience (TCN), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 7. Relationship of average noise and sensitivity of the stimulus response. 
%
% Created: Wed 10 Mar 2021, 10:32
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Sun 5 Jun 2022, 14:15
% Last edited by: Martin Esparza-Iaizzo

%% Add supplementary function and data path

addpath(genpath('./supple'));
addpath(genpath('./data_stimuli'));

%% Load data

% Load sigma = 0.02 data. UP
tmp = load('avg_up2_osc.mat');
sigma_low.avg_up_osc = tmp.avg_up_osc;
tmp = load('avg_up2_bis.mat');
sigma_low.avg_up_bis = tmp.avg_up_bis;
tmp = load('post_up2_bis.mat');
sigma_low.post_up_bis = tmp.post_up_bis;
tmp = load('post_up2_osc.mat');
sigma_low.post_up_osc = tmp.post_up_osc;

% Load sigma = 0.02 data. DOWN
tmp = load('avg_down2_osc.mat');
sigma_low.avg_down_osc = tmp.avg_down_osc;
tmp = load('avg_down2_bis.mat');
sigma_low.avg_down_bis = tmp.avg_down_bis;
tmp = load('post_down2_bis.mat');
sigma_low.post_down_bis = tmp.post_down_bis;
tmp = load('post_down2_osc.mat');
sigma_low.post_down_osc = tmp.post_down_osc;

% Load sigma = 0.25 data. UP
tmp = load('avg_up_osc.mat');
sigma_high.avg_up_osc = tmp.avg_up_osc;
tmp = load('avg_up_bis.mat');
sigma_high.avg_up_bis = tmp.avg_up_bis;
tmp = load('post_up_bis.mat');
sigma_high.post_up_bis = tmp.post_up_bis;
tmp = load('post_up_osc.mat');
sigma_high.post_up_osc = tmp.post_up_osc;

% Load sigma = 0.25 data. DOWN
tmp = load('avg_down_osc.mat');
sigma_high.avg_down_osc = tmp.avg_down_osc;
tmp = load('avg_down_bis.mat');
sigma_high.avg_down_bis = tmp.avg_down_bis;
tmp = load('post_down_bis.mat');
sigma_high.post_down_bis = tmp.post_down_bis;
tmp = load('post_down_osc.mat');
sigma_high.post_down_osc = tmp.post_down_osc;

%% Fitting. Sigma = 0.25

% Establish duration cut-off
dur = 10; 

% Clear UPs
[clean_osc_up,idx_osc_up] = clearCases(sigma_high.post_up_osc, size(sigma_high.post_up_osc) , dur);
sigma_high.avg_up_osc = sigma_high.avg_up_osc(idx_osc_up);
[clean_bis_up,idx_bis_up] = clearCases(sigma_high.post_up_bis, size(sigma_high.post_up_bis) , dur);
sigma_high.avg_up_bis = sigma_high.avg_up_bis(idx_bis_up);

% Clear Downs
[clean_osc_down,idx_osc_down] = clearCases(sigma_high.post_down_osc, size(sigma_high.post_down_osc) , dur);
sigma_high.avg_down_osc = sigma_high.avg_down_osc(idx_osc_down);
[clean_bis_down,idx_bis_down] = clearCases(sigma_high.post_down_bis, size(sigma_high.post_down_bis) , dur);
sigma_high.avg_down_bis = sigma_high.avg_down_bis(idx_bis_down);

% Generate p1, p2, p3, p4 and sensitivity values.
[sigma_high.p1_osc_up, sigma_high.p2_osc_up, sigma_high.p3_osc_up, sigma_high.p4_osc_up, sigma_high.sens_osc_up] = sigmoidRegression(clean_osc_up,'up');
[sigma_high.p1_bis_up, sigma_high.p2_bis_up, sigma_high.p3_bis_up, sigma_high.p4_bis_up, sigma_high.sens_bis_up] = sigmoidRegression(clean_bis_up,'up');
[sigma_high.p1_osc_down, sigma_high.p2_osc_down, sigma_high.p3_osc_down, sigma_high.p4_osc_down, sigma_high.sens_osc_down] = sigmoidRegression(clean_osc_down,'dw');
[sigma_high.p1_bis_down, sigma_high.p2_bis_down, sigma_high.p3_bis_down, sigma_high.p4_bis_down, sigma_high.sens_bis_down] = sigmoidRegression(clean_bis_down,'dw');

%% Fitting. Sigma = 0.02

% Establish duration cut-off
dur = 10; 

% Clear UPs
[clean_osc_up,idx_osc_up] = clearCases(sigma_low.post_up_osc, size(sigma_low.post_up_osc) , dur);
sigma_low.avg_up_osc = sigma_low.avg_up_osc(idx_osc_up);
[clean_bis_up,idx_bis_up] = clearCases(sigma_low.post_up_bis, size(sigma_high.post_up_bis) , dur);
sigma_low.avg_up_bis = sigma_low.avg_up_bis(idx_bis_up);

% Clear Downs
[clean_osc_down,idx_osc_down] = clearCases(sigma_low.post_down_osc, size(sigma_low.post_down_osc) , dur);
sigma_low.avg_down_osc = sigma_low.avg_down_osc(idx_osc_down);
[clean_bis_down,idx_bis_down] = clearCases(sigma_low.post_down_bis, size(sigma_low.post_down_bis) , dur);
sigma_low.avg_down_bis = sigma_low.avg_down_bis(idx_bis_down);

% Generate p1, p2, p3, p4 and sensitivity values.
[sigma_low.p1_osc_up, sigma_low.p2_osc_up, sigma_low.p3_osc_up, sigma_low.p4_osc_up, sigma_low.sens_osc_up] = sigmoidRegression(clean_osc_up,'up');
[sigma_low.p1_bis_up, sigma_low.p2_bis_up, sigma_low.p3_bis_up, sigma_low.p4_bis_up, sigma_low.sens_bis_up] = sigmoidRegression(clean_bis_up,'up');
[sigma_low.p1_osc_down, sigma_low.p2_osc_down, sigma_low.p3_osc_down, sigma_low.p4_osc_down, sigma_low.sens_osc_down] = sigmoidRegression(clean_osc_down,'dw');
[sigma_low.p1_bis_down, sigma_low.p2_bis_down, sigma_low.p3_bis_down, sigma_low.p4_bis_down, sigma_low.sens_bis_down] = sigmoidRegression(clean_bis_down,'dw');

%% Plotting
f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17 20]*1.75);
fontSize = 16;
MarkerSize = 12;
LineW = 1.0;

ax1 = axes('Position',[0.1 0.67 0.4 0.22]);
ax1.PositionConstraint = 'innerposition';
plot(sigma_high.avg_down_bis,abs(sigma_high.sens_bis_down),'o','LineWidth',LineW,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(sigma_high.avg_down_osc,abs(sigma_high.sens_osc_down),'o','LineWidth',LineW,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{DOWN}');
xlim([-0.8 0.8])
ylim([1 10e5])
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')


ax2 = axes('Position',[0.575 0.67 0.4 0.22]);
ax2.PositionConstraint = 'innerposition';
plot(sigma_high.avg_up_bis,abs(sigma_high.sens_bis_up),'o','LineWidth',LineW,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(sigma_high.avg_up_osc,abs(sigma_high.sens_osc_up),'o','LineWidth',LineW,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{UP}');
xlim([-0.8 0.8])
yticks([1 100 10000 1000000])
yticklabels({''})
ylim([1 10e5])
set(ax2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')


ax1 = axes('Position',[0.1 0.325 0.4 0.22]);
ax1.PositionConstraint = 'innerposition';
plot(sigma_low.avg_down_bis,abs(sigma_low.sens_bis_down),'o','LineWidth',LineW,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(sigma_low.avg_down_osc,abs(sigma_low.sens_osc_down),'o','LineWidth',LineW,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{DOWN}');
xlim([-0.08 0.08])
ylim([1 10e5])
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')


ax2 = axes('Position',[0.575 0.325 0.4 0.22]);
ax2.PositionConstraint = 'innerposition';
plot(sigma_low.avg_up_bis,abs(sigma_low.sens_bis_up),'o','LineWidth',LineW,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(sigma_low.avg_up_osc,abs(sigma_low.sens_osc_up),'o','LineWidth',LineW,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{UP}');
xlim([-0.08 0.08])
yticks([1 100 10000 1000000])
ylim([1 10e5])
yticklabels({''})
set(ax2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')

%% EXPORT

cd '/Volumes/GoogleDrive-101271366273470520077/My Drive/PaperBelen/Figures/Temp figures'
set(f,'Renderer','Painter')
exportgraphics(gcf,'fig7.pdf','Resolution',300,'BackgroundColor','none')

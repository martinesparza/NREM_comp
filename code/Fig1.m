%% Figures @NREM_comp
%
% Theoretical and Computational Neuroscience (TCN), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 1. Computational model and phase plane dynamics
%
% Created: Mon 16 Nov 2020, 10:05
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Sun 22 May 2022, 17:15
% Last edited by: Martin Esparza-Iaizzo
% 

%% Add supplementary function path

addpath(genpath('./supple'))

%% Generate nullclines

sample_r = 0:0.001:1; % Sample r vector to draw nullclines. 

% Define global variables from Levenstein et al., 2019
x0 = 5; 
r0 = 0.5;
k = 15;
b = 1;

% Oscillatory Regime
w = 6; I = 2.5;
nul_a_1 = 1./(1 + exp(-k.*(sample_r-r0))); % Nullcline 1
nul_r_1 = (1/b) * (w.*sample_r + I - x0 + log(((1./sample_r) - 1))); % Nullcline 2

% Bistable
w = 6.3; I = 2.35;
nul_a_2 = 1./(1 + exp(-k.*(sample_r-r0)));
nul_r_2 = (1/b) * (w.*sample_r + I - x0 + log(((1./sample_r) - 1)));

% Down
w = 6; I = 2.4;
nul_a_3 = 1./(1 + exp(-k.*(sample_r-r0)));
nul_r_3 = (1/b) * (w.*sample_r + I - x0 + log(((1./sample_r) - 1)));

% Up
w = 6.15; I = 2.5;
nul_a_4 = 1./(1 + exp(-k.*(sample_r-r0)));
nul_r_4 = (1/b) * (w.*sample_r + I - x0 + log(((1./sample_r) - 1)));

%% Plotting

fontSize = 16; 

f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 19.05 5]*1.75);

ax1_1 = axes('Position',[0.5 0.25 0.21 0.6]);
ax1_1.PositionConstraint = 'innerposition';
ax1_1.FontName = 'Arial';
plot(sample_r,nul_a_1,'Color',[.55 .55 .55],'LineWidth',3); hold on;
plot(sample_r,nul_r_1,'k','LineWidth',3);
ylim([0 1])
yticks([0 1]); yticklabels({'\color{black}0','\color{black}1'})
xticks([0 1]); xticklabels({'\color{black}0','\color{black}1'})
set(ax1_1,'FontSize',fontSize,'FontName','Arial','Box','on','LineWidth',2,'XColor','#FF44C8','YColor','#FF44C8')
xlabel('\color{black}\itr','FontName','Cambria')
ylabel('\color{black}\ita','FontName','Cambria')

ax1_2 = axes('Position',[0.74 0.25 0.21 0.6]);
ax1_2.PositionConstraint = 'innerposition';
plot(sample_r,nul_r_2,'k','LineWidth',3); hold on;
plot(sample_r,nul_a_2,'Color',[.55 .55 .55],'LineWidth',3); 
ylim([0 1])
yticks([]); yticklabels({''})
xticks([0 1]); xticklabels({'\color{black}0','\color{black}1'})
set(ax1_2,'FontSize',fontSize,'FontName','Arial','Box','on','LineWidth',2,'XColor','#027EDC','YColor','#027EDC')
xlabel('\color{black}\itr','FontName','Cambria')

%% Save figure ––––––– Uncomment and edit to save to personalised location 

% cd '/Users/martinesparzaiaizzo/Desktop/@belen/Figuras_LIMPIAS/METHODS/figures_temp/'
% set(f,'Renderer','Painter')
% exportgraphics(f,'fig1.pdf','Resolution',300,'BackgroundColor','none')


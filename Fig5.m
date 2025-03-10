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

%% New data load

% sigma high
cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/variables/sigma = 0.25'
files = dir('*.mat');
for i=1:length(files)
    load(files(i).name);
end

sigma_high.avg_up_osc = avg_up_osc;
sigma_high.avg_up_bis = avg_up_bis;
sigma_high.post_up_bis = post_up_bis(:,11);
sigma_high.post_up_osc = post_up_osc(:,11);
sigma_high.avg_down_osc = avg_down_osc;
sigma_high.avg_down_bis = avg_down_bis;
sigma_high.post_down_bis = post_down_bis(:,11);
sigma_high.post_down_osc = post_down_osc(:,11);

% sigma low
cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/variables/sigma = 0.0182'
files = dir('*.mat');
for i=1:length(files)
    load(files(i).name);
end

sigma_low.avg_up_osc = avg_up_osc;
sigma_low.avg_up_bis = avg_up_bis;
sigma_low.post_up_bis = post_up_bis;
sigma_low.post_up_osc = post_up_osc;
sigma_low.avg_down_osc = avg_down_osc;
sigma_low.avg_down_bis = avg_down_bis;
sigma_low.post_down_bis = post_down_bis;
sigma_low.post_down_osc = post_down_osc;
%% Load data

% Point to data location
cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/final_data_00182/UP'
% Load data
tmp = load('avg_up2_osc.mat');
sigma_low.avg_up_osc = tmp.avg_up_osc;
tmp = load('avg_up2_bis.mat');
sigma_low.avg_up_bis = tmp.avg_up_bis;
tmp = load('post_up2_bis.mat');
sigma_low.post_up_bis = tmp.post_up_bis;
tmp = load('post_up2_osc.mat');
sigma_low.post_up_osc = tmp.post_up_osc;


cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/final_data_00182/DOWN'
tmp = load('avg_down2_osc.mat');
sigma_low.avg_down_osc = tmp.avg_down_osc;
tmp = load('avg_down2_bis.mat');
sigma_low.avg_down_bis = tmp.avg_down_bis;
tmp = load('post_down2_bis.mat');
sigma_low.post_down_bis = tmp.post_down_bis;
tmp = load('post_down2_osc.mat');
sigma_low.post_down_osc = tmp.post_down_osc;

cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/sigma = 0.07/0.25/UP - 0.25'
tmp = load('avg_up_osc.mat');
sigma_high.avg_up_osc = tmp.avg_up_osc;
tmp = load('avg_up_bis.mat');
sigma_high.avg_up_bis = tmp.avg_up_bis;
tmp = load('post_up_bis.mat');
sigma_high.post_up_bis = tmp.post_up_bis;
tmp = load('post_up_osc.mat');
sigma_high.post_up_osc = tmp.post_up_osc;

cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/sigma = 0.07/0.25/DOWN - 0.25'
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
set(f,'units','centimeters','Position',[mp(2,1)+50 mp(2,1)+50 17 20]*1.75);
fontSize = 16;
MarkerSize = 10;

% A
ax1 = axes('Position',[0.1 0.67 0.32 0.25]);
ax1.PositionConstraint = 'innerposition';
plot(sigma_high.avg_down_bis,sigma_high.post_down_bis,'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on; 
plot(sigma_high.avg_down_osc,sigma_high.post_down_osc,'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8')
ylim([1 10000])
xlim([-0.8 0.8])
xlabel('<\xi>_t')
ylabel('Duration (AU)')
yticklabels({'10^0','','10^2','','10^4'})
yticks([1 10 100 1000 10000])
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
ax1_1 = axes('Position',[0.42 0.67 0.1 0.25]);
ax1_1.PositionConstraint = 'innerposition';
h = histogram(log10(sigma_high.post_down_bis), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#027EDC', 'EdgeAlpha', 0, 'FaceAlpha', 0.5); hold on
h.Orientation = 'horizontal';
h = histogram(log10(sigma_high.post_down_osc), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#FF44C8', 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
h.Orientation = 'horizontal';
ylim([0 4])
axis off

% B
ax2 = axes('Position',[0.56 0.67 0.32 0.25]);
ax2.PositionConstraint = 'innerposition';
plot(sigma_high.avg_up_bis,sigma_high.post_up_bis,'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on; 
plot(sigma_high.avg_up_osc,sigma_high.post_up_osc,'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8')
ylim([1 10000])
xlim([-0.8 0.8])
xlabel('<\xi>_t')
%ylabel('Up Duration (AU)')
yticklabels({'10^0','','10^2','','10^4'})
yticks([1 10 100 1000 10000])
set(ax2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
ax2_1 = axes('Position',[0.88 0.67 0.1 0.25]);
ax2_1.PositionConstraint = 'innerposition';
h = histogram(log10(sigma_high.post_up_bis), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#027EDC', 'EdgeAlpha', 0, 'FaceAlpha', 0.5); hold on
h.Orientation = 'horizontal';
h = histogram(log10(sigma_high.post_up_osc), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#FF44C8', 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
h.Orientation = 'horizontal';
ylim([0 4])
axis off

% C
ax3 = axes('Position',[0.1 0.3 0.32 0.25]);
ax3.PositionConstraint = 'innerposition';
plot(sigma_low.avg_down_bis,sigma_low.post_down_bis(:,11),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on; 
plot(sigma_low.avg_down_osc,sigma_low.post_down_osc(:,11),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8')
ylim([1 10000])
xlim([-0.08 0.08])
xlabel('<\xi>_t')
ylabel('Duration (AU)')
yticklabels({'10^0','','10^2','','10^4'})
yticks([1 10 100 1000 10000])
set(ax3,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
ax3_1 = axes('Position',[0.42 0.3 0.1 0.25]);
ax3_1.PositionConstraint = 'innerposition';
h = histogram(log10(sigma_low.post_down_bis(:,11)), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#027EDC', 'EdgeAlpha', 0, 'FaceAlpha', 0.5); hold on
h.Orientation = 'horizontal';
h = histogram(log10(sigma_low.post_down_osc(:,11)), 10, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#FF44C8', 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
h.Orientation = 'horizontal';
ylim([0 4])
axis off

% D
ax4 = axes('Position',[0.56 0.3 0.32 0.25]);
ax4.PositionConstraint = 'innerposition';
plot(sigma_low.avg_up_bis,sigma_low.post_up_bis(:,11),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on; 
plot(sigma_low.avg_up_osc,sigma_low.post_up_osc(:,11),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8')
ylim([1 10000])
xlim([-0.08 0.08])
xlabel('<\xi>_t')
yticklabels({'10^0','','10^2','','10^4'})
yticks([1 10 100 1000 10000])
set(ax4,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
ax4_1 = axes('Position',[0.88 0.3 0.1 0.25]);
ax4_1.PositionConstraint = 'innerposition';
h = histogram(log10(sigma_low.post_up_bis(:,11)), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#027EDC', 'EdgeAlpha', 0, 'FaceAlpha', 0.5); hold on
h.Orientation = 'horizontal';
h = histogram(log10(sigma_low.post_up_osc(:,11)), 10, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#FF44C8', 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
h.Orientation = 'horizontal';
ylim([0 4])
axis off

%% Save figure ––––––– Uncomment and edit to save to personalised location

% cd '/Users/martinesparzaiaizzo/Library/CloudStorage/GoogleDrive-martineladio.esparza01@alumni.upf.edu/My Drive/PaperBelen/Figures/Temp figures'
% set(f,'Renderer','Painter')
% exportgraphics(gcf,'fig5.pdf','Resolution',300,'BackgroundColor','none')


%% Stats –– Pearson correlation

% % Clear variables
% clear rho_sigma1_DOWN pval_sigma1_DOWN rho_sigma1_UP pval_sigma1_UP
% clear rho_sigma2_DOWN pval_sigma2_DOWN rho_sigma2_UP pval_sigma2_UP
% 
% 
% % Sigma 1 – DOWN
% [rho_sigma1_DOWN.bis,pval_sigma1_DOWN.bis] = corr(sigma_high.avg_down_bis,sigma_high.post_down_bis); % Bistable
% [rho_sigma1_DOWN.osc,pval_sigma1_DOWN.osc] = corr(sigma_high.avg_down_osc,sigma_high.post_down_osc); % Oscillatory
% 
% % Sigma 1 - UP
% [rho_sigma1_UP.bis,pval_sigma1_UP.bis] = corr(sigma_high.avg_up_bis,sigma_high.post_up_bis); % Bistable
% [rho_sigma1_UP.osc,pval_sigma1_UP.osc] = corr(sigma_high.avg_up_osc,sigma_high.post_up_osc); % Oscillatory
% 
% % Sigma 2 – DOWN
% [rho_sigma2_DOWN.bis,pval_sigma2_DOWN.bis] = corr(sigma_low.avg_down_bis,sigma_low.post_down_bis(:,11)); % Bistable
% [rho_sigma2_DOWN.osc,pval_sigma2_DOWN.osc] = corr(sigma_low.avg_down_osc,sigma_low.post_down_osc(:,11)); % Oscillatory
% 
% % Sigma 2 - UP
% [rho_sigma2_UP.bis,pval_sigma2_UP.bis] = corr(sigma_low.avg_up_bis,sigma_low.post_up_bis(:,11)); % Bistable
% [rho_sigma2_UP.osc,pval_sigma2_UP.osc] = corr(sigma_low.avg_up_osc,sigma_low.post_up_osc(:,11)); % Oscillatory


%% Levene tests
% sigma_high.p_down = vartestn([log(sigma_high.post_down_bis), log(sigma_high.post_down_osc)],'TestType','LeveneAbsolute');
% sigma_high.p_up = vartestn([log(sigma_high.post_up_bis), log(sigma_high.post_up_osc)],'TestType','LeveneAbsolute');
% 
% sigma_low.p_down = vartestn([log(sigma_low.post_down_bis(:,11)), log(sigma_low.post_down_osc(:,11))],'TestType','LeveneAbsolute');
% sigma_low.p_up = vartestn([log(sigma_low.post_up_bis(:,11)), log(sigma_low.post_up_osc(:,11))],'TestType','LeveneAbsolute');
% 

%% Anova

% Define variables
y = [sigma_high.post_up_bis; sigma_high.post_down_bis;...
    sigma_high.post_up_osc; sigma_high.post_down_osc;...
    sigma_low.post_up_bis(:,11); sigma_low.post_down_bis(:,11);...
    sigma_low.post_up_osc(:,11); sigma_low.post_down_osc(:,11)];

g1 = [repmat({'bis'}, 600,1); repmat({'osc'}, 600,1);...
    repmat({'bis'}, 600,1); repmat({'osc'}, 600,1)];

g2 = [repmat({'high noise'}, 1200,1); repmat({'low noise'}, 1200,1)];
g3 = [repmat({'up'}, 300,1); repmat({'down'}, 300,1);...
    repmat({'up'}, 300,1); repmat({'down'}, 300,1);...
    repmat({'up'}, 300,1); repmat({'down'}, 300,1);...
    repmat({'up'}, 300,1); repmat({'down'}, 300,1)];

% Run anova
[~, ~, stats] = anovan(y,{g1 g2},'model','interaction','varnames',{'g1','g2'});

% Run multcompare
[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);
tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A")=gnames(tbl.("Group A"));
tbl.("Group B")=gnames(tbl.("Group B"));

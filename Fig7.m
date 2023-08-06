%% Fig 4 Results

%% Sigma = 0.25
clear all
cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/variables/sigma = 0.25'
files = dir('*.mat');
for i=1:length(files)
    load(files(i).name);
end

cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/codeTower'

% clear all
dur = 10; 

% UPs
[clean_osc_up,idx_osc_up] = clearCases(post_up_osc, size(post_up_osc) , dur);
avg_up_osc = avg_up_osc(idx_osc_up);

[clean_bis_up,idx_bis_up] = clearCases(post_up_bis, size(post_up_bis) , dur);
avg_up_bis = avg_up_bis(idx_bis_up);

% Downs
[clean_osc_down,idx_osc_down] = clearCases(post_down_osc, size(post_down_osc) , dur);
avg_down_osc = avg_down_osc(idx_osc_down);

[clean_bis_down,idx_bis_down] = clearCases(post_down_bis, size(post_down_bis) , dur);
avg_down_bis = avg_down_bis(idx_bis_down);


% Generate p1, p2, p3, p4 and sensitivity values.
[p1_osc_up, p2_osc_up, p3_osc_up, p4_osc_up, sens_osc_up] = sigmoidRegression(clean_osc_up,'up');
[p1_bis_up, p2_bis_up, p3_bis_up, p4_bis_up, sens_bis_up] = sigmoidRegression(clean_bis_up,'up');
[p1_osc_down, p2_osc_down, p3_osc_down, p4_osc_down, sens_osc_down] = sigmoidRegression(clean_osc_down,'dw');
[p1_bis_down, p2_bis_down, p3_bis_down, p4_bis_down, sens_bis_down] = sigmoidRegression(clean_bis_down,'dw');

%% Plotting
f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 17 20]*1.75);
fontSize = 16;
MarkerSize = 10;

ax1 = axes('Position',[0.1 0.7 0.375 0.25]);
ax1.PositionConstraint = 'innerposition';
plot(avg_down_bis,abs(sens_bis_down),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(avg_down_osc,abs(sens_osc_down),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{DOWN}');
xlim([-0.8 0.8])
ylim([1 10e5])
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
ax1_1 = axes('Position', get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax1_1.YAxis.Visible = 'off';
ax1_1.XLim = [-1.6 1.6];
ax1_1.XTick = [-1 0 1];
set(ax1_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')


ax2 = axes('Position',[0.575 0.7 0.375 0.25]);
ax2.PositionConstraint = 'innerposition';
plot(avg_up_bis,abs(sens_bis_up),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(avg_up_osc,abs(sens_osc_up),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{UP}');
xlim([-0.8 0.8])
ylim([1 10e5])
set(ax2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
ax1_1 = axes('Position', get(ax2,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax1_1.YAxis.Visible = 'off';
ax1_1.XLim = [-1.6 1.6];
ax1_1.XTick = [-1 0 1];
set(ax1_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')

tmp = NaN(size(sens_osc_down,1),2);
tmp(1:length(sens_bis_down),1) = log(abs(sens_bis_down)); tmp(1:length(sens_osc_down),2) = log(abs(sens_osc_down)); 
sigma_high.p_down = vartestn(tmp,'TestType','LeveneAbsolute');

tmp = NaN(size(sens_osc_up,1),2);
tmp(1:length(sens_bis_up),1) = log(abs(sens_bis_up)); tmp(1:length(sens_osc_up),2) = log(abs(sens_osc_up)); 
sigma_high.p_up = vartestn(tmp,'TestType','LeveneAbsolute');


%% Sigma = 0.02
cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/variables/sigma = 0.0182'
files = dir('*.mat');
for i=1:length(files)
    load(files(i).name);
end

cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/codeTower'

% clear all
dur = 10; 

% UPs
[clean_osc_up,idx_osc_up] = clearCases(post_up_osc, size(post_up_osc) , dur);
avg_up_osc = avg_up_osc(idx_osc_up);

[clean_bis_up,idx_bis_up] = clearCases(post_up_bis, size(post_up_bis) , dur);
avg_up_bis = avg_up_bis(idx_bis_up);

% Downs
[clean_osc_down,idx_osc_down] = clearCases(post_down_osc, size(post_down_osc) , dur);
avg_down_osc = avg_down_osc(idx_osc_down);

[clean_bis_down,idx_bis_down] = clearCases(post_down_bis, size(post_down_bis) , dur);
avg_down_bis = avg_down_bis(idx_bis_down);


% Generate p1, p2, p3, p4 and sensitivity values.
[p1_osc_up, p2_osc_up, p3_osc_up, p4_osc_up, sens_osc_up] = sigmoidRegression(clean_osc_up,'up');
[p1_bis_up, p2_bis_up, p3_bis_up, p4_bis_up, sens_bis_up] = sigmoidRegression(clean_bis_up,'up');
[p1_osc_down, p2_osc_down, p3_osc_down, p4_osc_down, sens_osc_down] = sigmoidRegression(clean_osc_down,'dw');
[p1_bis_down, p2_bis_down, p3_bis_down, p4_bis_down, sens_bis_down] = sigmoidRegression(clean_bis_down,'dw');

%% Plotting

ax1 = axes('Position',[0.1 0.27 0.375 0.25]);
ax1.PositionConstraint = 'innerposition';
plot(avg_down_bis,abs(sens_bis_down),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(avg_down_osc,abs(sens_osc_down),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{DOWN}');
xlim([-0.08 0.08])
ylim([1 10e5])
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
ax1_1 = axes('Position', get(ax1,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax1_1.YAxis.Visible = 'off';
ax1_1.XLim = [-.16 .16];
ax1_1.XTick = [-.1 0 .1];
set(ax1_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')


ax2 = axes('Position',[0.575 0.27 0.375 0.25]);
ax2.PositionConstraint = 'innerposition';
plot(avg_up_bis,abs(sens_bis_up),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(avg_up_osc,abs(sens_osc_up),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{UP}');
xlim([-0.08 0.08])
ylim([1 10e5])
set(ax2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
ax1_1 = axes('Position', get(ax2,'Position'), ...
    'XAxisLocation','top', ...
    'Color','none', ...
    'XColor','k');
ax1_1.YAxis.Visible = 'off';
ax1_1.XLim = [-.16 .16];
ax1_1.XTick = [-.1 0 .1];
set(ax1_1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')

sigma_low.p_down = vartestn([log(abs(sens_bis_down)), log(abs(sens_osc_down))],'TestType','LeveneAbsolute');
sigma_low.p_up = vartestn([log(abs(sens_bis_up)), log(abs(sens_osc_up))],'TestType','LeveneAbsolute');


% 
% 
% ax3 = axes('Position',[0.1 0.1 0.175 0.175]);
% ax3.PositionConstraint = 'innerposition';
% plot(clean_bis_down(:,11),p4_bis_down,'o','LineWidth',1.0,'MarkerSize',5,'Color','#027EDC'); hold on;
% plot(clean_osc_down(:,11),p4_osc_down,'o','LineWidth',1.0,'MarkerSize',5,'Color','#FF44C8');
% plot(sort(clean_bis_down(:,11)),sort(clean_bis_down(:,11)),'--k','LineWidth',1.5);
% xlim([0 1000])
% ylim([0 350])
% xlabel('Duration (AU)'); ylabel('p4_{DOWN}');
% set(ax3,'FontSize',9,'Box','on','LineWidth',1.5,'FontName','Arial','xscale','log','XMinorTick','off')
% 
%%





%% EXPORT
% 
% cd '/Users/martinesparzaiaizzo/Library/CloudStorage/GoogleDrive-martineladio.esparza01@alumni.upf.edu/My Drive/PaperBelen/Figures/Temp figures'
% exportgraphics(gcf,'fig7.pdf','Resolution',300,'BackgroundColor','none')
% 
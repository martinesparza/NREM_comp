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
avg_up_osc_high = avg_up_osc(idx_osc_up);

[clean_bis_up,idx_bis_up] = clearCases(post_up_bis, size(post_up_bis) , dur);
avg_up_bis_high = avg_up_bis(idx_bis_up);

% Downs
[clean_osc_down,idx_osc_down] = clearCases(post_down_osc, size(post_down_osc) , dur);
avg_down_osc_high = avg_down_osc(idx_osc_down);

[clean_bis_down,idx_bis_down] = clearCases(post_down_bis, size(post_down_bis) , dur);
avg_down_bis_high = avg_down_bis(idx_bis_down);


% Generate p1, p2, p3, p4 and sensitivity values.
[p1_osc_up, p2_osc_up, p3_osc_up, p4_osc_up, sens_osc_up_high] = sigmoidRegression(clean_osc_up,'up');
[p1_bis_up, p2_bis_up, p3_bis_up, p4_bis_up, sens_bis_up_high] = sigmoidRegression(clean_bis_up,'up');
[p1_osc_down, p2_osc_down, p3_osc_down, p4_osc_down, sens_osc_down_high] = sigmoidRegression(clean_osc_down,'dw');
[p1_bis_down, p2_bis_down, p3_bis_down, p4_bis_down, sens_bis_down_high] = sigmoidRegression(clean_bis_down,'dw');

%% Plotting
f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,1)+50 17 20]*1.75);
fontSize = 16;
MarkerSize = 10;

ax1 = axes('Position',[0.1 0.67 0.32 0.25]);
ax1.PositionConstraint = 'innerposition';
plot(avg_down_bis_high,abs(sens_bis_down_high),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(avg_down_osc_high,abs(sens_osc_down_high),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('$S$', 'Interpreter', 'latex');
xlim([-0.8 0.8])
ylim([1 10e5])
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
ax1_1 = axes('Position',[0.42 0.67 0.1 0.25]);
ax1_1.PositionConstraint = 'innerposition';
h = histogram(log10(abs(sens_bis_down_high)), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#027EDC', 'EdgeAlpha', 0, 'FaceAlpha', 0.5); hold on
h.Orientation = 'horizontal';
h = histogram(log10(abs(sens_osc_down_high)), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#FF44C8', 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
h.Orientation = 'horizontal';
ylim([0 6])
axis off

ax2 = axes('Position',[0.56 0.67 0.32 0.25]);
ax2.PositionConstraint = 'innerposition';
plot(avg_up_bis_high,abs(sens_bis_up_high),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(avg_up_osc_high,abs(sens_osc_up_high),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t');
xlim([-0.8 0.8])
ylim([1 10e5])
set(ax2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
% yticklabels({''})
ax2_1 = axes('Position',[0.88 0.67 0.1 0.25]);
ax2_1.PositionConstraint = 'innerposition';
h = histogram(log10(abs(sens_bis_up_high)), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#027EDC', 'EdgeAlpha', 0, 'FaceAlpha', 0.5); hold on
h.Orientation = 'horizontal';
h = histogram(log10(abs(sens_osc_up_high)), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#FF44C8', 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
h.Orientation = 'horizontal';
ylim([0 6])
axis off


%% Sigma = 0.02
cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/variables/sigma = 0.0182 corrected'
files = dir('*.mat');
for i=1:length(files)
    load(files(i).name);
end

cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/codeTower'

% clear all
dur = 10; 

% UPs
[clean_osc_up,idx_osc_up] = clearCases(post_up_osc, size(post_up_osc) , dur);
avg_up_osc_low = avg_up_osc(idx_osc_up);

[clean_bis_up,idx_bis_up] = clearCases(post_up_bis, size(post_up_bis) , dur);
avg_up_bis_low = avg_up_bis(idx_bis_up);

% Downs
[clean_osc_down,idx_osc_down] = clearCases(post_down_osc, size(post_down_osc) , dur);
avg_down_osc_low = avg_down_osc(idx_osc_down);

[clean_bis_down,idx_bis_down] = clearCases(post_down_bis, size(post_down_bis) , dur);
avg_down_bis_low = avg_down_bis(idx_bis_down);


% Generate p1, p2, p3, p4 and sensitivity values.
[p1_osc_up, p2_osc_up, p3_osc_up, p4_osc_up, sens_osc_up_low] = sigmoidRegression(clean_osc_up,'up');
[p1_bis_up, p2_bis_up, p3_bis_up, p4_bis_up, sens_bis_up_low] = sigmoidRegression(clean_bis_up,'up');
[p1_osc_down, p2_osc_down, p3_osc_down, p4_osc_down, sens_osc_down_low] = sigmoidRegression(clean_osc_down,'dw');
[p1_bis_down, p2_bis_down, p3_bis_down, p4_bis_down, sens_bis_down_low] = sigmoidRegression(clean_bis_down,'dw');

%% Plotting
ax1 = axes('Position',[0.1 0.3 0.32 0.25]);
ax1.PositionConstraint = 'innerposition';
plot(avg_down_bis_low,abs(sens_bis_down_low),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(avg_down_osc_low,abs(sens_osc_down_low),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('$S$', 'Interpreter', 'latex');
xlim([-0.08 0.08])
ylim([1 10e5])
set(ax1,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
ax1_1 = axes('Position',[0.42 0.3 0.1 0.25]);
ax1_1.PositionConstraint = 'innerposition';
h = histogram(log10(abs(sens_bis_down_low)), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#027EDC', 'EdgeAlpha', 0, 'FaceAlpha', 0.5); hold on
h.Orientation = 'horizontal';
h = histogram(log10(abs(sens_osc_down_low)), 10, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#FF44C8', 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
h.Orientation = 'horizontal';
ylim([0 6])
axis off

ax2 = axes('Position',[0.56 0.3 0.32 0.25]);
ax2.PositionConstraint = 'innerposition';
plot(avg_up_bis_low,abs(sens_bis_up_low),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#027EDC'); hold on;
plot(avg_up_osc_low,abs(sens_osc_up_low),'o','LineWidth',1.0,'MarkerSize',MarkerSize,'Color','#FF44C8');
xlabel('<\xi>_t');
xlim([-0.08 0.08])
ylim([1 10e5])
set(ax2,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')
% yticklabels({''})
ax2_1 = axes('Position',[0.88 0.3 0.1 0.25]);
ax2_1.PositionConstraint = 'innerposition';
h = histogram(log10(abs(sens_bis_up_low)), 15, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#027EDC', 'EdgeAlpha', 0, 'FaceAlpha', 0.5); hold on
h.Orientation = 'horizontal';
h = histogram(log10(abs(sens_osc_up_low)), 10, 'Normalization', 'count', 'DisplayStyle', 'bar',...
    'FaceColor', '#FF44C8', 'EdgeAlpha', 0, 'FaceAlpha', 0.5);
h.Orientation = 'horizontal';
ylim([0 6])
axis off

% 
% sigma_low.p_down = vartestn([log(abs(sens_bis_down)), log(abs(sens_osc_down))],'TestType','LeveneAbsolute');
% sigma_low.p_up = vartestn([log(abs(sens_bis_up)), log(abs(sens_osc_up))],'TestType','LeveneAbsolute');


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



%% EXPORT
% 
cd '/Users/martinesparzaiaizzo/Library/CloudStorage/GoogleDrive-martineladio.esparza01@alumni.upf.edu/My Drive/PaperBelen/Figures/Temp figures'
exportgraphics(gcf,'sup_sens_duration.pdf','Resolution',300,'BackgroundColor','none')


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
set(f,'units','centimeters','Position',[0 0 19.05 19]);

ax1 = axes('Position',[0.1 0.7 0.375 0.25]);
ax1.PositionConstraint = 'innerposition';
plot(avg_down_bis,abs(sens_bis_down),'o','LineWidth',1.0,'MarkerSize',5,'Color','#027EDC'); hold on;
plot(avg_down_osc,abs(sens_osc_down),'o','LineWidth',1.0,'MarkerSize',5,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{DOWN}');
xlim([-0.8 0.8])
ylim([1 10e5])
set(ax1,'FontSize',9,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')


ax2 = axes('Position',[0.575 0.7 0.375 0.25]);
ax2.PositionConstraint = 'innerposition';
plot(avg_up_bis,abs(sens_bis_up),'o','LineWidth',1.0,'MarkerSize',5,'Color','#027EDC'); hold on;
plot(avg_up_osc,abs(sens_osc_up),'o','LineWidth',1.0,'MarkerSize',5,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{UP}');
xlim([-0.8 0.8])
ylim([1 10e5])
set(ax2,'FontSize',9,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')




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
plot(avg_down_bis,abs(sens_bis_down),'o','LineWidth',1.0,'MarkerSize',5,'Color','#027EDC'); hold on;
plot(avg_down_osc,abs(sens_osc_down),'o','LineWidth',1.0,'MarkerSize',5,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{DOWN}');
xlim([-0.08 0.08])
ylim([1 10e5])
set(ax1,'FontSize',9,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')


ax2 = axes('Position',[0.575 0.27 0.375 0.25]);
ax2.PositionConstraint = 'innerposition';
plot(avg_up_bis,abs(sens_bis_up),'o','LineWidth',1.0,'MarkerSize',5,'Color','#027EDC'); hold on;
plot(avg_up_osc,abs(sens_osc_up),'o','LineWidth',1.0,'MarkerSize',5,'Color','#FF44C8');
xlabel('<\xi>_t'); ylabel('Sensitivity_{UP}');
xlim([-0.08 0.08])
ylim([1 10e5])
set(ax2,'FontSize',9,'Box','on','LineWidth',1.5,'FontName','Arial','yscale','log','YMinorTick','off')

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
% ax4 = axes('Position',[0.3 0.1 0.175 0.175]);
% ax4.PositionConstraint = 'innerposition';
% figure
% plot(avg_up_bis,p2_bis_up,'o','LineWidth',1.0,'MarkerSize',5,'Color','#027EDC'); hold on;
% plot(avg_up_osc,p2_osc_up,'o','LineWidth',1.0,'MarkerSize',5,'Color','#FF44C8');
% plot(sort(clean_bis_down(:,11)),sort(clean_bis_down(:,11)),'--k','LineWidth',1.5);
% xlim([0 1000])
% ylim([0 350])
% yticks([]); yticklabels({''})
% xlabel('Noise'); ylabel('p4_{UP}');
% set(ax4,'FontSize',9,'Box','on','LineWidth',1.5,'FontName','Arial','XMinorTick','off')

%% EXPORT
% 
% cd '/Users/martinesparzaiaizzo/Desktop/@belen/Figuras_LIMPIAS/RESULTS'
% exportgraphics(gcf,'fig4.eps','Resolution',300,'BackgroundColor','none')


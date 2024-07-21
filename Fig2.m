%% Figures @NREM_comp
%
% Theoretical and Computational Neuroscience (TCN), Centre for Brain and Cognition, 
% Universitat Pompeu Fabra, Barcelona, Spain. 
%
% Figure 2. Extended model and stimulus. 
%
% Created: Fri 22 Jan 2021, 14:31
% Author: Martin Esparza-Iaizzo
% 
% Last edited: Sun 22 May 2022, 17:19
% Last edited by: Martin Esparza-Iaizzo
%
%% Plotting

duration_s = 10; % Stimulus duration
dt = 0.01; % Time interval 

fontSize = 16; 

f = figure; 
mp = get(0, 'MonitorPositions');
set(f,'units','centimeters','Position',[mp(end,1)+50 mp(end,2)+50 19.05 20]*1.75);

% Fig. 2B. Example of Poisson Spike Generation. 
ax = axes('Position',[0.5 0.8 0.2 0.09]);
ax.PositionConstraint = 'innerposition';

freqs = 5; % Frequency of spiking
time = 5; % Duration of Poisson Spike train
trials = 3; % Number of sample trials 
[spikeMat, tVec] = poissonSpikeGen(freqs, time, trials, dt);

plotRaster(spikeMat, tVec,1,'#DAA520');
xlim([-1 6])
yticks([1 3]); yticklabels({'1','3'})
xticks([0 5]); xticklabels({'0','5'})
xlabel('Time (AU)')
ylabel('Trials')
ax.YAxis.TickLength = [0 0];
set(ax,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')

ax0 = axes('Position',[0.5 0.9 0.2 0.05]);
ax0.PositionConstraint = 'innerposition';
plot([zeros(1,100) ones(1,500) zeros(1,100)], 'LineWidth',1.5,'Color','#DAA520')
xlim([0 700])
ylim([-0.25 1.25])
xticks([100 600]); xticklabels({''})
yticks([0 1]); yticklabels({'0','5'})
ylabel('\lambda')
ax0.YAxis.TickLength = [0 0];
set(ax0,'FontSize',fontSize,'Box','off','LineWidth',1.5,'FontName','Arial')

% Fig. 2C. Alpha-function. 
[~,~,t_tmp] = generateStimuli(duration_s,dt,1,freqs);
alfa = 1/0.1;
alfa_function = alfa^2 .* t_tmp(1:250) .* exp(-alfa.*t_tmp(1:250));

ax3 = axes('Position',[0.775 0.8 0.175 0.15]);
ax3.PositionConstraint = 'innerposition';
plot(t_tmp(1:250),alfa_function,'LineWidth',1.5,'Color','k');
xlim([0 0.75])
yticks([0 max(alfa_function)]); yticklabels({'0','1'})
xticks([0 0.75]); xticklabels({'0','3'})
xlabel('Time (AU)');
ylabel('{\alpha - function}')
set(ax3,'FontSize',fontSize,'Box','on','LineWidth',1.5,'FontName','Arial')


% Fig. 2D. Spike train + alpha-function convolution
% Spiking frequency = 5
freqs = [5];
[spike_train,spikeMat,tVec] = generateStimuli(duration_s,dt,1,freqs);

ax1_1 = axes('Position',[0.1 0.54 0.85 0.075]);
ax1_1.PositionConstraint = 'innerposition';
plot(tVec,spike_train(2,:),'Color','#DAA520','LineWidth',1.5); hold on;
ylim([0 1]);
yticks([0 1]); yticklabels({'0','1'});
xticks([]); xticklabels({''})
xlabel('');
set(ax1_1,'FontSize',fontSize,'Box','off','LineWidth',1.5,'FontName','Arial')
ylabel('\its(t)','FontName','Arial')

ax1_2 = axes('Position',[0.1 0.615 0.85 0.07]);
ax1_2.PositionConstraint = 'innerposition';
ax1_2.XAxis.Visible = 'off';
plotRaster(spikeMat, tVec,1.5,'#DAA520');
yticks([]); yticklabels({''})
xticks([]); xticklabels({''})
ylim([0.25 2])
set(ax1_2,'FontSize',fontSize,'Box','off','LineWidth',1.5,'FontName','Arial')


% Spiking frequency = 10
freqs = [10];
[spike_train,spikeMat,tVec] = generateStimuli(duration_s,dt,1,freqs);

ax2_1 = axes('Position',[0.1 0.36 0.85 0.075]);
ax2_1.PositionConstraint = 'innerposition';
plot(tVec,spike_train(2,:),'Color','#DAA520','LineWidth',1.5); hold on;
ylim([0 1]);
yticks([0 1]); yticklabels({'0','1'});
xticks([0 10]); xticklabels({'0','10'})
xlabel('');
xlabel('Time (AU)')
set(ax2_1,'FontSize',fontSize,'Box','off','LineWidth',1.5,'FontName','Arial')
ylabel('\its(t)','FontName','Arial')

ax2_2 = axes('Position',[0.1 0.435 0.85 0.07]);
ax2_2.PositionConstraint = 'innerposition';
ax2_2.XAxis.Visible = 'off';
plotRaster(spikeMat, tVec,1.5,'#DAA520');
yticks([]); yticklabels({''})
xticks([]); xticklabels({''})
ylim([0.25 2])
set(ax2_2,'FontSize',fontSize,'Box','off','LineWidth',1.5,'FontName','Arial')

% Spiking frequency = 50
freqs = [50];
[spike_train,spikeMat,tVec] = generateStimuli(duration_s,dt,1,freqs);

ax2_1 = axes('Position',[0.1 0.1 0.85 0.075]);
ax2_1.PositionConstraint = 'innerposition';
plot(tVec,spike_train(2,:),'Color','#DAA520','LineWidth',1.5); hold on;
ylim([0 1]);
yticks([0 1]); yticklabels({'0','1'});
xticks([0 10]); xticklabels({'0','10'})
xlabel('');
xlabel('Time (AU)')
set(ax2_1,'FontSize',fontSize,'Box','off','LineWidth',1.5,'FontName','Arial')
ylabel('\its(t)','FontName','Arial')

ax2_2 = axes('Position',[0.1 0.175 0.85 0.07]);
ax2_2.PositionConstraint = 'innerposition';
ax2_2.XAxis.Visible = 'off';
plotRaster(spikeMat, tVec,1.5,'#DAA520');
yticks([]); yticklabels({''})
xticks([]); xticklabels({''})

ylim([0.25 2])
set(ax2_2,'FontSize',fontSize,'Box','off','LineWidth',1.5,'FontName','Arial')


%% Save figure ––––––– Uncomment and edit to save to personalised location
% 

str = sprintf('/Users/martinesparzaiaizzo/Library/CloudStorage/GoogleDrive-martineladio.esparza01@alumni.upf.edu/My Drive/PaperBelen/Figures/Temp figures/fig2.pdf');
set(f,'Renderer','Painter')
exportgraphics(gcf,str,'Resolution',300,'BackgroundColor','none')

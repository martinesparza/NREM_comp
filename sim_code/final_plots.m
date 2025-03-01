%% Final plots

% cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/variables/sigma = 0.25'
cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/variables/sigma = 0.0182'
files = dir('*.mat');
for i=1:length(files)
    load(files(i).name);
end



cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/codeTower'

% clear all

[clean_out_o,clean_idx_o,fractions_o, idxs ] = clearCases(post_down_osc, size(post_up_osc) ,10);
avg_down_osc = avg_down_osc(clean_idx_o);

[clean_out_b,clean_idx_b,fractions_b, idxs ] = clearCases(post_down_bis, size(post_up_bis) ,10);
avg_down_bis = avg_down_bis(clean_idx_b);


%% sensitivity

size_inp = size(clean_out_o);

p3_values_o =zeros(size_inp(1),1); sens_o = zeros(size_inp(1),1); idx_o = [];p4_values_o = zeros(size_inp(1),1);
p1_values_o = zeros(size_inp(1),1); p2_values_o = zeros(size_inp(1),1);

f2 = -100:10:100;
for i = 1:size_inp(1)
    
    y = clean_out_o(i,:);
    fun = @(p,x) p(4) + (p(1)-p(4))./(1+exp(-p(2)*(x-p(3)))); 
    x = linspace(-100,100);
    d = f2; %frequencies
    p0 = zeros(1,4);
    p0(1) = max([y(1) y(end)]);
    
    p0(4) = min([y(1) y(end)]);
    slope = diff(y)/5;
    [~, loc] = max(abs(slope));
    p0(2) = slope(loc);
    p0(3) = (d(loc)+d(loc+1))/2;

    lb = [0 -inf -100 0];
    ub = [Inf Inf 100 Inf];
    
    error_fun = @(p_new)(p_new(4) + (p_new(1)-p_new(4))./(1+exp(-p_new(2)*(d-p_new(3)))))-y;

    p_fitted = lsqnonlin(error_fun,p0,lb,ub);
    
    if (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4)) < 0
        p3_values_o(i) = p_fitted(3);
        p4_values_o(i) = p_fitted(4);
        p1_values_o(i) = p_fitted(1);
        p2_values_o(i) = p_fitted(2);
        sens_o(i) = (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4));
    end
    if (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4)) >= 0
        p3_values_o(i) = NaN;
        p4_values_o(i) = NaN;
        p1_values_o(i) = NaN;
        p2_values_o(i) = NaN;
        sens_o(i) = NaN;
    end

        
%     figure;
%     plot(d,y,'ko',x,fun(p_fitted,x));
%     legend('Data','Best fit')
%     xlabel('Frequency')
%     xlim([-100 100])
%     ylabel('Duration (s)')
%     str = sprintf('Sensitivity %f',p2_values_b(i));
%     title(str);
%     

end

size_inp = size(clean_out_b);

p3_values_b =zeros(size_inp(1),1); sens_b = zeros(size_inp(1),1); idx_b = [];p4_values_b = zeros(size_inp(1),1);
p1_values_b = zeros(size_inp(1),1); p2_values_b = zeros(size_inp(1),1);

f2 = -100:10:100;
for i = 1:size_inp(1)
    
    y = clean_out_b(i,:);
    fun = @(p,x) p(4) + (p(1)-p(4))./(1+exp(-p(2)*(x-p(3)))); 
    x = linspace(-100,100);
    d = f2; %frequencies
    p0 = zeros(1,4);
    p0(1) = max([y(1) y(end)]);
    
    p0(4) = min([y(1) y(end)]);
    slope = diff(y)/5;
    [~, loc] = max(abs(slope));
    p0(2) = slope(loc);
    p0(3) = (d(loc)+d(loc+1))/2;

    lb = [0 -inf -100 0];
    ub = [Inf Inf 100 Inf];
    
    error_fun = @(p_new)(p_new(4) + (p_new(1)-p_new(4))./(1+exp(-p_new(2)*(d-p_new(3)))))-y;

    p_fitted = lsqnonlin(error_fun,p0,lb,ub);
    
    if (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4)) < 0
        p3_values_b(i) = p_fitted(3);
        p4_values_b(i) = p_fitted(4);
        p1_values_b(i) = p_fitted(1);
        p2_values_b(i) = p_fitted(2);
        sens_b(i) = (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4));
    end
    if (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4)) >= 0
        p3_values_b(i) = NaN;
        p4_values_b(i) = NaN;
        p1_values_b(i) = NaN;
        p2_values_b(i) = NaN;
        sens_b(i) = NaN;
    end

        

end

%%
% P3 
% cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/Marvin/figs/sigma = 0.25/up'
cd '/Users/martinesparzaiaizzo/Desktop/UPF 19-20/2o trimestre/Practicas CBC/Marvin/figs/sigma = 0.0182/up'

figure; 
plot(avg_down_bis,p3_values_b,'o','LineWidth',2.0,'MarkerSize',12); hold on;
plot(avg_down_osc,p3_values_o,'o','LineWidth',2.0,'MarkerSize',12);
xlabel('<\xi>_t'); ylabel('p3');
% title('p3 vs Noise');
% xlim([-0.05 0.05])
ylim([-100 100])
legend('Bistable','Oscillatory')
legend('Bistable','Oscillatory')
set(gca,'FontSize',24,'box','off','LineWidth',2)
set(gcf,'Position',[0 0 800 600])
% saveas(gcf,'p3.png')

%P1
figure; 
plot(avg_down_bis,p1_values_b,'o','LineWidth',2.0,'MarkerSize',12); hold on;
plot(avg_down_osc,p1_values_o,'o','LineWidth',2.0,'MarkerSize',12);
xlabel('<\xi>_t'); ylabel('p1');
% title('p1 vs Noise');
% xlim([-0.05 0.05])
ylim([1 1000])
legend('Bistable','Oscillatory')
set(gca,'FontSize',24,'box','off','LineWidth',2)
set(gca,'yscale','log')
set(gcf,'Position',[0 0 800 600])
% saveas(gcf,'p1.png')

% P2
figure; 
plot(avg_down_bis,p2_values_b,'o','LineWidth',2.0,'MarkerSize',12); hold on;
plot(avg_down_osc,p2_values_o,'o','LineWidth',2.0,'MarkerSize',12);
xlabel('<\xi>_t'); ylabel('p2');
% xlim([-0.05 0.05])
ylim([0.1 1000])
legend('Bistable','Oscillatory')
set(gca,'FontSize',24,'box','off','LineWidth',2)
set(gcf,'Position',[0 0 800 600])
set(gca,'yscale','log')
% saveas(gcf,'p2.png')

% P4
figure; 
plot(avg_down_bis,p4_values_b,'o','LineWidth',2.0,'MarkerSize',12); hold on;
plot(avg_down_osc,p4_values_o,'o','LineWidth',2.0,'MarkerSize',12);
xlabel('<\xi>_t'); ylabel('p4');
% title('p4 vs Noise');
% xlim([-0.05 0.05])
ylim([1 100])
legend('Bistable','Oscillatory')
set(gca,'FontSize',24,'box','off','LineWidth',2)
set(gca,'yscale','log')
set(gcf,'Position',[0 0 800 600])
% saveas(gcf,'p4.png')

% Sens
figure; 
plot(avg_down_bis,log(sens_b),'o','LineWidth',2.0,'MarkerSize',12); hold on;
plot(avg_down_osc,log(sens_o),'o','LineWidth',2.0,'MarkerSize',12);
% plot(linspace(0,10,10000),spike_trains_10sec([1:2:19],:),'LineWidth',2.0)
xlabel('<\xi>_t'); ylabel('log(slope)');
% title('Slope vs. Noise.');
set(gca,'FontSize',24,'box','off','LineWidth',2)
set(gcf,'Position',[0 0 800 600])
legend('Bistable','Oscillatory')
% xlim([-0.05 0.05])
ylim([-2 14])
% saveas(gcf,'slope.png')

% Dur
figure; 
plot(avg_down_bis,post_up_bis(clean_idx_b,11),'o','LineWidth',2.0,'MarkerSize',12); hold on;
plot(avg_down_osc,post_up_osc(clean_idx_o,11),'o','LineWidth',2.0,'MarkerSize',12);
xlabel('<\xi>_t'); ylabel('Durations');
% title('Duration vs. Noise');
% xlim([-0.05 0.05])
ylim([1 10000])
set(gca,'yscale','log')
legend('Bistable','Oscillatory')
set(gca,'FontSize',24,'box','off','LineWidth',2)
set(gcf,'Position',[0 0 800 600])
% saveas(gcf,'dur.png')



%% Parallel Plot

% bistable_stimuli = [avg_stimuli_bis p3_values_b log(abs(sens_b))];
% bistable_up = [avg_down_bis p3_values_b log(abs(sens_b))];
% 
% osc_stimuli = [avg_stimuli_osc p3_values_o log(abs(sens_o))];
% osc_up = [avg_down_osc p3_values_o log(abs(sens_o))];
% 
% stimuli = [bistable_stimuli; osc_stimuli];
% up = [bistable_up; osc_up];
% 
% labels_stimuli = cell(length(bistable_stimuli) + length(osc_stimuli),1);
% labels_stimuli(1:length(bistable_stimuli),1) = {'Bistable'};
% labels_stimuli(length(bistable_stimuli)+1:end,1) = {'Oscillatory'};
% 
% labels_up = cell(length(bistable_up) + length(osc_stimuli),1);
% labels_up(1:length(bistable_up),1) = {'Bistable'};
% labels_up(length(bistable_up)+1:end,1) = {'Oscillatory'};
% 
% figure;
% coorddata = [1:3];
% p = parallelplot(stimuli, 'GroupData',labels_stimuli, 'CoordinateData',coorddata);
% title('Average over stimulus')
% set(gca,'FontSize',24,'LineWidth',2)
% p.CoordinateTickLabels = {'Noise','p3','Sensitivity'};
% % set(gcf,'Position',[0 0 800 600])
% 
% figure;
% p = parallelplot(up, 'GroupData',labels_up, 'CoordinateData',coorddata);
% title('Average over DOWN')
% set(gca,'FontSize',24,'LineWidth',2)
% p.CoordinateTickLabels = {'Noise','p3','Sensitivity'};




%% Mas figuras. p1 p2 p3 p4 vs duration

% [p3_b,idx3] = rmmissing(p3_values_b);
% [p4_b,idx4] = rmmissing(p4_values_b);
% [p1_b,idx1] = rmmissing(p1_values_b);
% [p2_b,idx2] = rmmissing(p2_values_b);
% [sb,idxb] = rmmissing(sens_b);
% 
% [p3_o,idx3] = rmmissing(p3_values_o);
% [p4_o,idx4] = rmmissing(p4_values_o);
% [p1_o,idx1] = rmmissing(p1_values_o);
% [p2_o,idx2] = rmmissing(p2_values_o);
% [so,idxo] = rmmissing(sens_o);
% 
% 
% figure; 
% plot(clean_out_b((find(idxb == 0)),11),sb,'o','LineWidth',2.0,'MarkerSize',12); hold on;
% plot(clean_out_o((find(idxo == 0)),11),so,'o','LineWidth',2.0,'MarkerSize',12); hold on;
% xlabel('Durations'); ylabel('Slope');
% title('DOWN.');
% % xlim([-0.6 0.8])
% set(gca,'xscale','log')
% legend('Bistable','Oscillatory')
% set(gca,'FontSize',24,'box','off','LineWidth',2)
% set(gcf,'Position',[0 0 800 600])
% % 
% 
% figure;  
% plot(clean_out_b((find(idxb == 0)),11),p1_b,'o','LineWidth',2.0,'MarkerSize',12); hold on;
% plot(clean_out_o((find(idxo == 0)),11),p1_o,'o','LineWidth',2.0,'MarkerSize',12); hold on;
% plot(sort(clean_out_b(:,11)),sort(clean_out_b(:,11)),'--k','LineWidth',1.5);
% xlabel('Durations'); ylabel('p1');
% title('DOWN.');
% % xlim([-0.6 0.8])
% set(gca,'xscale','log')
% legend('Bistable','Oscillatory','Identity')
% set(gca,'FontSize',24,'box','off','LineWidth',2)
% set(gcf,'Position',[0 0 800 600])
% 
% figure; 
% plot(clean_out_b((find(idxb == 0)),11),p2_b,'o','LineWidth',2.0,'MarkerSize',12); hold on;
% plot(clean_out_o((find(idxo == 0)),11),p2_o,'o','LineWidth',2.0,'MarkerSize',12); hold on;
% xlabel('Durations'); ylabel('p2');
% title('DOWN.');
% % xlim([-0.6 0.8])
% set(gca,'xscale','log')
% legend('Bistable','Oscillatory')
% set(gca,'FontSize',24,'box','off','LineWidth',2)
% set(gcf,'Position',[0 0 800 600])
% 
% figure; 
% plot(clean_out_b((find(idxb == 0)),11),p3_b,'o','LineWidth',2.0,'MarkerSize',12); hold on;
% plot(clean_out_o((find(idxo == 0)),11),p3_o,'o','LineWidth',2.0,'MarkerSize',12); hold on;
% xlabel('Durations'); ylabel('p3');
% title('DOWN.');
% % xlim([-0.6 0.8])
% set(gca,'xscale','log')
% legend('Bistable','Oscillatory')
% set(gca,'FontSize',24,'box','off','LineWidth',2)
% set(gcf,'Position',[0 0 800 600])
% 
% 
% figure; 
% plot(clean_out_b((find(idxb == 0)),11),p4_b,'o','LineWidth',2.0,'MarkerSize',12); hold on;
% plot(clean_out_o((find(idxo == 0)),11),p4_o,'o','LineWidth',2.0,'MarkerSize',12); hold on;
% plot(sort(clean_out_b(:,11)),sort(clean_out_b(:,11)),'--k','LineWidth',1.5);
% xlabel('Durations'); ylabel('p4');
% title('DOWN.');
% % xlim([-0.6 0.8])
% set(gca,'xscale','log')
% legend('Bistable','Oscillatory','Identity')
% set(gca,'FontSize',24,'box','off','LineWidth',2)
% set(gcf,'Position',[0 0 800 600])


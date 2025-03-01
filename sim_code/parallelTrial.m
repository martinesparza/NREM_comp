%% Definitions

clear all
clc
close all
% r = mean firing rate
% a = adaptation

% Arbitrary units for the time constants in order to avoid dimensions
tau_r = 1;
tau_a = 25; 

% Fixed values for sigmoidal function
x0 = 5; 
r0 = 0.5;
k = 15;

% Weight values can change for different regimes
% Values for oscillatory regime
b = 1;
w = 6.3;
I = 2.35;

% Noise function is employed with Ornstein - Uhlenbeck method

duration = 1000; %s
save_dt = 0.01; % Precision on saving values in noise function and in solving the ode
dt = 0.01; % must be the same as save_dt
theta =0.05; %About comparable...
sigma = 0.25;
numsignals = 1;

[noise, noise_t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);

% New input. Stimuli

% amplitude = 1;
% duration_s = 1000; %10s 
% start = randi(100000 - duration);
% interval = start:(start + duration);
% 
% stimuli = zeros(1,100001);
% stimuli(60000:70000) = amplitude; %Amplitude of stimul


%%

tic
n = 200;
A = 500;
a = zeros(n);
parfor i = 1:n
    a(i) = max(abs(eig(rand(A))));
end
toc

%% Initialize cell memory UP

start = 60000;
stimuli = zeros(1,100001);
max_sim = 100; % total number of simulation (100) 
amplitude_vec = 0.5:0.1:2.5;
duration_vec = 50:10:50; % This means that duration will be simulated from 10s to 100s in intervals of 10s
memory_cell_up = initializePar(duration_vec,amplitude_vec,max_sim,dt);

noise_memory = zeros(max_sim,length(memory_cell_up),length(amplitude_vec));
ic_memory = zeros(max_sim,2,length(amplitude_vec));
 % counter for every row.


% While loop for temporal average. It can be for "r" or for "a".

% UP CASE
tic

for j = 1:length(duration_vec)
duration_s = duration_vec(j)/dt;
margin = 100;
margin_i = margin/dt;


    for i = 1:length(amplitude_vec)
        stimuli(start:start + duration_s) = amplitude_vec(i);

        parfor parall = 1:max_sim
            n_up = 0;
            while n_up < 1
                [noise, noise_t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
                ic = rand(2,1);

                [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0,stimuli), [0:dt:duration], ic);
                r_s = ra_s(:,1);
                %a_s = ra_s(:,2);

                [rs_thresh,rs_cross,rs_bihist,rs_diptest] = BimodalThresh(r_s,'Schimdt');
                %[a_thresh,a_cross,a_bihist,a_diptest] = BimodalThresh(a_s,'Schmidt');

                dim_up_r = length(rs_cross.upints);
                N = histcounts(start,reshape(rs_cross.upints',[],1)); % Vamos a comprobar si el inicio esta en un intervalo up

                if mod(find(N == 1),2) ~= 0
                    n_up = 1
                    memory_cell_up(parall,:,i) = r_s((start-margin_i):(start + duration_s +margin_i));
                    noise_memory(parall,:,i) = noise((start-margin_i):(start + duration_s +margin_i));
                    ic_memory(parall,:,i) = ic'
                end
            end
        end
    end
end
toc

% aux_t = -margin:(duration_s+margin); %Auxiliary temporal vector 
% figure;
% plot(aux_t,memory_up(1:n_up,:)');


%% Initialize cell memory DOWN

start = 60000;
stimuli = zeros(1,100001);
max_sim = 100; % total number of simulation (100) 
amplitude_vec = 0.5:0.1:2.5;
duration_vec = 80:10:80;
0; % This means that duration will be simulated from 10s to 100s in intervals of 10s
memory_cell_down = initializePar(duration_vec,amplitude_vec,max_sim,dt);

noise_memory = zeros(max_sim,length(memory_cell_down),length(amplitude_vec));
ic_memory = zeros(max_sim,2,length(amplitude_vec));
% While loop for temporal average. It can be for "r" or for "a".

% DOWN CASE
tic

for j = 1:length(duration_vec)
duration_s = duration_vec(j)/dt;
margin = 100;
margin_i = margin/dt;


    for i = 1:length(amplitude_vec)
        stimuli(start:start + duration_s) = amplitude_vec(i);

        parfor parall = 1:max_sim
            n_down = 0;
            while n_down < 1
                [noise, noise_t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
                ic = rand(2,1);

                [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0,stimuli), [0:dt:duration], ic);
                r_s = ra_s(:,1);
                %a_s = ra_s(:,2);

                [rs_thresh,rs_cross,rs_bihist,rs_diptest] = BimodalThresh(r_s,'Schimdt');
                %[a_thresh,a_cross,a_bihist,a_diptest] = BimodalThresh(a_s,'Schmidt');

                dim_up_r = length(rs_cross.downints);
                N = histcounts(start,reshape(rs_cross.downints',[],1)); % Vamos a comprobar si el inicio esta en un intervalo up

                if mod(find(N == 1),2) ~= 0
                    n_down = 1;
                    memory_cell_down(parall,:,i) = r_s((start-margin_i):(start + duration_s +margin_i));
                    noise_memory(parall,:,i) = noise((start-margin_i):(start + duration_s +margin_i));
                    ic_memory(parall,:,i) = ic';
                end
            end
        end
    end
end
toc

% aux_t = -margin:(duration_s+margin); %Auxiliary temporal vector 
% figure;
% plot(aux_t,memory_up(1:n_up,:)');




%% PDFS

tic
sim = 2000;
pdfs = zeros(sim,duration/dt + 1);
parfor i = 1:sim
    [noise, noise_t] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
    ic = rand(2,1);
    [t,ra] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0), [0:dt:duration], ic);
    pdfs(i,:) = ra(:,1);
    
end


vec = zeros((duration/dt+1)*sim,1);
vec(1:100001,1) = pdfs(1,:);
for i = 1:(sim-1)
    vec((((duration/dt+1)*i)+1):((duration/dt+1)*i)+100001,1) = pdfs(i+1,:);
end

[vec_thresh,vec_cross,vec_bihist,vec_diptest] = BimodalThresh(vec,'Schimdt');

% UP

size_up = length(vec_cross.upints);
pdf_dur_up = zeros(1,size_up);

parfor i=1:size_up
    pdf_dur_up(1,i) = vec_cross.upints(i,2) -  vec_cross.upints(i,1);
    
end
pdf_dur_up = pdf_dur_up*dt;



% DOWN

size_down = length(vec_cross.downints);
pdf_dur_down = zeros(1,size_down);

parfor i=1:size_down
    pdf_dur_down(1,i) = vec_cross.downints(i,2) -  vec_cross.downints(i,1);
    
end
pdf_dur_down = pdf_dur_down*dt;
[histo_down,bins_down] = hist(pdf_dur_down);

figure
plot(bins_down,histo_down/sum(histo_down),'-o')
title('Stimulus duration 100s')
xlabel('Time (s)')
ylabel('Probability')


toc
%% Recreate


[t,ra_rec] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0,stimuli), [0:dt:duration], ic(1,:,1));
r_rec = ra_s(:,1);
a_rec = ra_s(:,2);

%% Figures

x = linspace(0,(25001-1)*dt,25001) - 100;
figure;
sgtitle('Stimulation over UP intervals')
subplot(2,1,1)
plot(x,memory_cell_up(1,:,1),'LineWidth',4); hold on;
plot(x,noise_memory(1,:,1),'LineWidth',0.01)
legend('Mean firing rate','Noise')
xlabel('Time (AU)')
ylabel('Amplitude')
subplot(2,1,2)
plot(x,memory_cell_up(27,:,1),'LineWidth',4); hold on;
plot(x,noise_memory(27,:,1),'LineWidth',0.01)
legend('Mean firing rate','Noise')
xlabel('Time (AU)')
ylabel('Amplitude')

%% Durations

% UP
tic
thresh = 0.5499;
pre_post_duration = zeros(100,2,21);

for b = 1:21
    for a = 1:100
        vec = memory_cell_down(a,:,b);
        for i = (100/dt + 1):-1:1
            if vec(i) < thresh
                break 
            end
        end


        for n = (100/dt + 1):length(vec)
            if vec(n) < thresh
                break 
            end
        end

        pre_post_duration(a,1,b) = 100/dt + 1 - i;
        pre_post_duration(a,2,b) = n - (100/dt + 1);
    end
end
toc


figure
pre_post_duration = pre_post_duration*dt;
for i = [1 10 21]
    [histo,bins] = hist(pre_post_duration(:,1,i) + pre_post_duration(:,2,i));
    plot(bins,histo/sum(histo),'-o')
    hold on
end

% % %  DOWN

tic
thresh = 0.5499;
pre_post_ten = zeros(100,2,21);

for b = 1:21
    for a = 1:100
        vec = memory_ten(a,:,b);
        for i = (100/dt + 1):-1:1
            if vec(i) > thresh
                break 
            end
        end


        for n = (100/dt + 1):length(vec)
            if vec(n) > thresh
                break 
            end
        end

        pre_post_ten(a,1,b) = 100/dt + 1 - i;
        pre_post_ten(a,2,b) = n - (100/dt + 1);
    end
end
toc
pre_post_ten = pre_post_ten*dt;

% figure

% for i = [1 10 21]
%     [histo,bins] = hist(pre_post_duration(:,1,i) + pre_post_duration(:,2,i));
%     plot(bins,histo/sum(histo),'-o')
%     hold on
% end
% title('Stimulus duration 100s')
% xlabel('Time (s)')
% ylabel('Probability')
% legend('Amplitude = 0.5', 'Amplitude = 1.5', 'Amplitude = 2.5')

figure
for i = [1 10 21]
    plot(pre_post_ten(:,1,i),pre_post_ten(:,2,i),'o');
    hold on
end  
title('Pre vs Post. 50s Stimulation')
xlabel('pre duration)')
ylabel('post duration')
legend('Amplitude = 0.5', 'Amplitude = 1.5', 'Amplitude = 2.5')


figure
sgtitle('Pre vs Post. 50s Stimulation')
subplot(3,1,1)
plot(pre_post_ten(:,1,1),pre_post_ten(:,2,1),'o');
title('Amplitude 0.5')
subplot(3,1,2)
plot(pre_post_ten(:,1,10),pre_post_ten(:,2,10),'o');
title('Amplitude 1.5')
subplot(3,1,3)
plot(pre_post_ten(:,1,21),pre_post_ten(:,2,21),'o');
title('Amplitude 2.5')


%% New plots

% Noise vs. post duration. 50s
noise_average_fifty = zeros(100,3);
for i= 1:100
    noise_average_fifty(i,1) = mean(noise_fifty(i,10000:15000,1));
    noise_average_fifty(i,2) = mean(noise_fifty(i,10000:15000,11));
    noise_average_fifty(i,3) = mean(noise_fifty(i,10000:15000,21));
end
figure
plot(noise_average_fifty(:,1),pre_post_fifty(:,2,1),'o'); hold on;
plot(noise_average_fifty(:,2),pre_post_fifty(:,2,11),'o'); hold on;
plot(noise_average_fifty(:,3),pre_post_fifty(:,2,21),'o'); 
title('Noise average vs Post. 50s Stimulation')
xlabel('Noise average')
ylabel('post duration')
legend('Amplitude = 0.5', 'Amplitude = 1.5', 'Amplitude = 2.5')

figure
sgtitle('Pre vs Post. 50s Stimulation')
subplot(3,1,1)
plot(noise_average_fifty(:,1),pre_post_fifty(:,2,1),'o');
title('Amplitude 0.5')
subplot(3,1,2)
plot(noise_average_fifty(:,2),pre_post_fifty(:,2,11),'o');
title('Amplitude 1.5')
subplot(3,1,3)
plot(noise_average_fifty(:,3),pre_post_fifty(:,2,21),'o');
title('Amplitude 2.5')


% Noise vs. post duration. 10s
noise_average_ten = zeros(100,3);
for i= 1:100
    noise_average_ten(i,1) = mean(noise_ten(i,10000:11000,1));
    noise_average_ten(i,2) = mean(noise_ten(i,10000:11000,11));
    noise_average_ten(i,3) = mean(noise_ten(i,10000:11000,21));
end
figure
plot(noise_average_ten(:,1),pre_post_ten(:,2,1),'o'); hold on;
plot(noise_average_ten(:,2),pre_post_ten(:,2,11),'o'); hold on;
plot(noise_average_ten(:,3),pre_post_ten(:,2,21),'o'); 
title('Noise average vs Post. 10s Stimulation')
xlabel('Noise average')
ylabel('post duration')
legend('Amplitude = 0.5', 'Amplitude = 1.5', 'Amplitude = 2.5')

figure
sgtitle('Pre vs Post. 10s Stimulation')
subplot(3,1,1)
plot(noise_average_ten(:,1),pre_post_ten(:,2,1),'o');
title('Amplitude 0.5')
subplot(3,1,2)
plot(noise_average_ten(:,2),pre_post_ten(:,2,11),'o');
title('Amplitude 1.5')
subplot(3,1,3)
plot(noise_average_ten(:,3),pre_post_ten(:,2,21),'o');
title('Amplitude 2.5')

% Noise vs. post duration. 100s
noise_average_hundred = zeros(100,3);
for i= 1:100
    noise_average_hundred(i,1) = mean(noise_hundred(i,10000:20000,1));
    noise_average_hundred(i,2) = mean(noise_hundred(i,10000:20000,11));
    noise_average_hundred(i,3) = mean(noise_hundred(i,10000:20000,21));
end
figure
plot(noise_average_hundred(:,1),pre_post_hundred(:,2,1),'o'); hold on;
plot(noise_average_hundred(:,2),pre_post_hundred(:,2,11),'o'); hold on;
plot(noise_average_hundred(:,3),pre_post_hundred(:,2,21),'o'); 
title('Noise average vs Post. 100s Stimulation')
xlabel('Noise average')
ylabel('post duration')
legend('Amplitude = 0.5', 'Amplitude = 1.5', 'Amplitude = 2.5')

figure
sgtitle('Pre vs Post. 100s Stimulation')
subplot(3,1,1)
plot(noise_average_hundred(:,1),pre_post_hundred(:,2,1),'o');
title('Amplitude 0.5')
subplot(3,1,2)
plot(noise_average_hundred(:,2),pre_post_hundred(:,2,11),'o');
title('Amplitude 1.5')
subplot(3,1,3)
plot(noise_average_hundred(:,3),pre_post_hundred(:,2,21),'o');
title('Amplitude 2.5')


%% Ventanas de ruido. Valores fijos

% ii_cc = rand(2,1);
% [noise_test, noise_t_test,aux,randnums] = OUNoise(theta,sigma,duration,dt,save_dt,numsignals);
% trial_noise = zeros(100,100001);
% for i = 1:100
%     trial_noise(i,:) = noise_test;
% end

% options = odeset('OutputFcn',@(t,y,flag) myOutputFcn(t,y,flag,w,b,I,noise_test,noise_t,x0));
[t,ra] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise_test,noise_t,x0,k,r0), [0:dt:duration], ii_cc);%,options);
r = ra(:,1);
plot(t,noise_test,'LineWidth',1);
hold on; plot(t, new_memory(:,1:3,1),'LineWidth',0.01)
xlim([300 550])
legend('Firing rate','noise')
%plot(t,r_inf(1:(end-1)),'-o');


%%

x = linspace(0,(21001-1)*dt,21001) - 100;
start = 44000;
duration_s = 45/dt;
stimuli = zeros(1,100001);
% window_duration = duration - start*0.01;

% final_avg = zeros(3,100);
final_post = zeros(100,1);
j = 1;
% new_memory = zeros(100001,100,3);

% for start = 34000:6000:46000
    
    stimuli = zeros(1,100001);
    new_noise = noise_test(1:start); 
    complete_window = zeros(100001,1);
    pre = zeros(100,1); 
    avg = zeros(1,100); 
    post = zeros(100,1);
    
    stimuli(start:start + duration_s) = 1;
    window_duration = duration - start*0.01;
    
    
    for i = 1:100
%         [window,window_t] = OUNoiseWindow(theta,0.5,window_duration,dt,save_dt,numsignals,noise_test(start+1));
%         complete_window = [new_noise; window];
%         new_memory(:,i,j) = complete_window;
        
%         [t,ra] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,new_memory(:,i,j),noise_t,x0,k,r0), [0:dt:duration], ii_cc);
%         r = ra(:,1)';
%         [pre,post_no_stimuli] = durationMeasurement(r',start);
        
        [t,ra_s] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,new_memory(:,i,j),noise_t,x0,k,r0,stimuli), [0:dt:duration], ii_cc);
        r_s = ra_s(:,1)';
        [pre,post_stimuli] = durationMeasurement(r_s',start);
        
        post(i,1) = post_stimuli; 
%         avg(i) = mean(complete_window(start:start+post_no_stimuli));
    end
%     post = post*dt;
%     final_avg(j,:) = avg;
    final_post(:,j) = post;
    j = j + 1;
    
% end


%%

%new_new = [final_post(:,1)/180 final_post(:,2)/120 final_post(:,3)/60];
% figure
% for i = 1:3
%     plot(final_avg(i,:),new_new(:,i),'o'); hold on;
%     title('Noise vs. Post. 10s. 0.5 Amplitude')
%     xlabel('Noise average')
%     ylabel('Post duration (s)')
%     legend('520s onset','550s onset','580s onset')
% end



for scenario = 1
%     scenario = 3;

% aux = log(final_post(:,scenario));
% [n,xbins,ybins,xedges,y_log_edges_two] = hist2d(final_avg(scenario,:),aux);
%     [n,xbins,ybins,xedges,y_linear_edges] = hist2d(final_avg(scenario,:),final_post(:,scenario));
    
    aux = final_post(:,scenario) .* nostimuli(:,scenario);
    
    figure
    plot(final_avg(scenario,:),aux,'^','LineWidth',1,'MarkerFaceColor',[1 .5 1],'Color',[1 .5 1],'MarkerSize',20 ); hold on;
    plot(final_avg(scenario,:),nostimuli(:,scenario),'^','LineWidth',1,'MarkerEdgeColor',[1 .5 1],'MarkerSize',20)
%     histogram2(final_avg(scenario,:)',final_post(:,scenario),overall_x,exp(overall_edges),'FaceColor','flat','ShowEmptyBins','on','Normalization','probability');
%     hist3([final_avg(scenario,:)' final_post(:,scenario)],'Edges',{overall_x exp(overall_edges)},'CDataMode','auto','FaceColor','interp');
        
%         if scenario == 1
%             title('Noise vs. Post. No stimuli. 10s Onset')
%         end
%         if scenario == 2
%             title('Noise vs. Post. No stimuli. 70s Onset')
%         end
%         if scenario == 3
%             title('Noise vs. Post. No stimuli. 130s Onset')
%         end
    xlabel('noise')
    ylabel('Duration (s)')
    set(gca, 'YScale', 'log')%,'YMinorTick','off','YMinorGrid','off')
    legend('stimuli','no stimuli')
    set(gcf,'Position',[500 500 700 400])
    set(gca,'FontSize',22,'Box','off','LineWidth',1.5)
    xlim([-0.8 0.6])
%     colorbar
%     yticks(test)
%     xticks(round(overall_x,2))
%     view([0 90])
    %ylim([0 4])
%     if scenario == 1
%         saveas(gcf,'nostimuli_10s.tif')  
%     end
%     if scenario == 2
%         saveas(gcf,'nostimuli_70s.tif') 
%     end
%     if scenario == 3
%         saveas(gcf,'nostimuli_130s.tif') 
%     end
    
end


%% More figures. Down. 

i = [-0.8 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2];
j = [-0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.6];
cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/test/15-05/Down/NewData/data'
load('final_scaled5.mat')
load('final_scaled10.mat')
load('final_scaled25.mat')
load('final_scaled50.mat')
cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/test/15-05/Down'
load ('overall_avg.mat')
cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/test/15-05/Down/NewData/Figures10s'


for n = 1:8 
    

ten = final_avg(1,:);
seventy = final_avg(2,:);
ind_ten = find (ten >= i(n) & ten < j(n)); %length(ind_ten)
ind_seventy = find (seventy >= 0 & seventy < 0.1);
mean_vector = zeros(4,4);

mean_vector(1,1) = mean(final_post_scaled5(ind_ten,1));
mean_vector(2,1) = std(final_post_scaled5(ind_ten,1));
mean_vector(3,1) = mean(final_post_scaled5(ind_seventy,2));
mean_vector(4,1) = std(final_post_scaled5(ind_seventy,2));


mean_vector(1,2) = mean(final_post_scaled10(ind_ten,1));
mean_vector(2,2) = std(final_post_scaled10(ind_ten,1));
mean_vector(3,2) = mean(final_post_scaled10(ind_seventy,2));
mean_vector(4,2) = std(final_post_scaled10(ind_seventy,2));

mean_vector(1,3) = mean(final_post_scaled25(ind_ten,1));
mean_vector(2,3) = std(final_post_scaled25(ind_ten,1));
mean_vector(3,3) = mean(final_post_scaled25(ind_seventy,2));
mean_vector(4,3) = std(final_post_scaled25(ind_seventy,2));


mean_vector(1,4) = mean(final_post_scaled50(ind_ten,1));
mean_vector(2,4) = std(final_post_scaled50(ind_ten,1));
mean_vector(3,4) = mean(final_post_scaled50(ind_seventy,2));
mean_vector(4,4) = std(final_post_scaled50(ind_seventy,2));


figure
%plot([5 10 25 50],mean_vector(1,:),'o','MarkerSize',12,'LineWidth',1,'MarkerFaceColor',[.12 0.56 1],'Color',[.12 0.56 1],'MarkerEdgeColor','k'); hold on;
errorbar([5 10 25 50],mean_vector(1,:),mean_vector(2,:),'o','MarkerSize',12,'LineWidth',1,'MarkerFaceColor',[.12 0.56 1],'Color',[.12 0.56 1],'MarkerEdgeColor','k')
%plot([5 10 25 50],mean_vector(3,:),'o','MarkerSize',12,'LineWidth',1,'MarkerFaceColor',[1 .5 1],'Color',[1 .5 1],'MarkerEdgeColor','k');


str = sprintf('Down. Average duration vs. Stimulus. Noise: [%g, %g)',i(n),j(n));
title(str)
xlabel('Area under the stimuli curve')
ylabel('Average duration (s)')
xlim([0 55])
ylim([0 200])
%legend('10s onset', '70s onset')
% sva = sprintf('10s_Noise_%g_%g.tif',i(n),j(n));
% saveas(gcf,sva)  

end

%% More figures. Up. 

i = [-0.8 -0.2 -0.1 0.05 0.15 0.2 0.25 0.35];
j = [-0.2 -0.1 0.05 0.15 0.2 0.25 0.35 0.6];
cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/test/29-05/UP/NewData/data'
load('final_scaled_5.mat')
load('final_scaled_10.mat')
load('final_scaled_25.mat')
load('final_scaled_50.mat')
cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/test/15-05/UP'
load ('overall_avg.mat')
% cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/test/29-05/UP/NewData/Figures10s'


for n = 1:8 
    

ten = final_avg(1,:);
seventy = final_avg(2,:);
ind_ten = find (ten >= i(n) & ten < j(n)); %length(ind_ten)
ind_seventy = find (seventy >= 0 & seventy < 0.1);
mean_vector = zeros(4,4);

mean_vector(1,1) = mean(final_scaled_5(ind_ten,1));
mean_vector(2,1) = std(final_scaled_5(ind_ten,1));
mean_vector(3,1) = mean(final_scaled_5(ind_seventy,2));
mean_vector(4,1) = std(final_scaled_5(ind_seventy,2));


mean_vector(1,2) = mean(final_scaled_10(ind_ten,1));
mean_vector(2,2) = std(final_scaled_10(ind_ten,1));
mean_vector(3,2) = mean(final_scaled_10(ind_seventy,2));
mean_vector(4,2) = std(final_scaled_10(ind_seventy,2));

mean_vector(1,3) = mean(final_scaled_25(ind_ten,1));
mean_vector(2,3) = std(final_scaled_25(ind_ten,1));
mean_vector(3,3) = mean(final_scaled_25(ind_seventy,2));
mean_vector(4,3) = std(final_scaled_25(ind_seventy,2));


mean_vector(1,4) = mean(final_scaled_50(ind_ten,1));
mean_vector(2,4) = std(final_scaled_50(ind_ten,1));
mean_vector(3,4) = mean(final_scaled_50(ind_seventy,2));
mean_vector(4,4) = std(final_scaled_50(ind_seventy,2));


figure
%plot([5 10 25 50],mean_vector(1,:),'o','MarkerSize',12,'LineWidth',1,'MarkerFaceColor',[.12 0.56 1],'Color',[.12 0.56 1],'MarkerEdgeColor','k'); hold on;
errorbar([5 10 25 50],mean_vector(1,:),mean_vector(2,:),'o','MarkerSize',12,'LineWidth',1,'MarkerFaceColor',[.12 0.56 1],'Color',[.12 0.56 1],'MarkerEdgeColor','k')
%plot([5 10 25 50],mean_vector(3,:),'o','MarkerSize',12,'LineWidth',1,'MarkerFaceColor',[1 .5 1],'Color',[1 .5 1],'MarkerEdgeColor','k');


str = sprintf('Down. Average duration vs. Stimulus. Noise: [%g, %g)',i(n),j(n));
title(str)
xlabel('Area under the stimuli curve')
ylabel('Average duration (s)')
xlim([0 55])
ylim([0 200])
%legend('10s onset', '70s onset')
% sva = sprintf('10s_Noise_%g_%g.tif',i(n),j(n));
% saveas(gcf,sva)  

end

%% Final figures?
clear all
clc

% UP amplitud fija
% cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/test/29-05/UP/Amplitud fija'
% load('final_scaled05.mat'); load('final_scaled15.mat'); load('final_scaled25.mat'); load('final_scaled35.mat'); load('final_scaled45.mat');
% load('final_scaled10.mat'); load('final_scaled20.mat'); load('final_scaled30.mat'); load('final_scaled40.mat'); load('final_scaled50.mat');

% UP duracion fija
% cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/test/29-05/UP/Duracion fija'
% load('final_post01.mat'); load('final_post02.mat'); load('final_post03.mat'); load('final_post04.mat'); load('final_post05.mat');
% load('final_post06.mat'); load('final_post07.mat'); load('final_post08.mat'); load('final_post09.mat'); load('final_post1.mat');

% UP
% cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/test/29-05/UP/Final Data'
% load('Up_amplitudesVariables.mat')
% load('Up_duracionesVariables.mat')
% load('avg_up.mat')
% load('nostimuli_up.mat')

% cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/test/29-05/Down/Final data'
% load('Down_amplitudesVariables.mat')
% load('Down_duracionesVariables.mat')
% load('avg_down.mat')
% load('nostimuli_down.mat')

cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/200s_Oscillatory'
load('Down_Oscillatory_amplitudeVariable.mat')
load('Down_Oscillatory_durationVariable.mat')
load('Down_Oscillatory_nostimuli.mat')
load('Up_Oscillatory_amplitudeVariable.mat')
load('Up_Oscillatory_durationVariable.mat')
load('Up_Oscillatory_nostimuli.mat')
load('oscillatory_avg_up_200.mat')
load('oscillatory_avg_down_200.mat')


i_up = [-0.8 -0.2 -0.1 0.05 0.15 0.2 0.25 0.35];
j_up = [-0.2 -0.1 0.05 0.15 0.2 0.25 0.35 0.6];

i_down = [-0.8 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2];
j_down = [-0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.6];

i_test = [-0.8 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4];
j_test = [-0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.6];

cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/Paper/Figures/Fig9/Oscillatory'

for n = 1
    mean_vector = zeros(8,11); % Fila de datos y fila de error
    ind_up = find (avg_up_200(:,1) >= i_test(n) & avg_up_200(:,1) < j_test(n));
    ind_down = find (avg_down_200(:,1) >= i_test(n) & avg_down_200(:,1) < j_test(n));

    mean_vector(1,1) = mean(Up_Oscillatory_nostimuli(ind_up));
    mean_vector(2,1) = std(Up_Oscillatory_nostimuli(ind_up))/sqrt(length(ind_up));
    mean_vector(3,1) = mean_vector(1,1);
    mean_vector(4,1) = mean_vector(2,1);

    mean_vector(5,1) = mean(Down_Oscillatory_nostimuli(ind_down));
    mean_vector(6,1) = std(Down_Oscillatory_nostimuli(ind_down))/sqrt(length(ind_down));
    mean_vector(7,1) = mean_vector(5,1);
    mean_vector(8,1) = mean_vector(6,1);

    for k = 2:11
        mean_vector(1,k) = mean(Up_Oscillatory_durationVariable(ind_up,k-1));
        mean_vector(2,k) = std(Up_Oscillatory_durationVariable(ind_up,k-1))/sqrt(length(ind_up));

        mean_vector(3,k) = mean(Up_Oscillatory_amplitudeVariable(ind_up,k-1));
        mean_vector(4,k) = std(Up_Oscillatory_amplitudeVariable(ind_up,k-1))/sqrt(length(ind_up));

        mean_vector(5,k) = mean(Down_Oscillatory_durationVariable(ind_down,k-1));
        mean_vector(6,k) = std(Down_Oscillatory_durationVariable(ind_down,k-1))/sqrt(length(ind_down));

        mean_vector(7,k) = mean(Down_Oscillatory_amplitudeVariable(ind_down,k-1));
        mean_vector(8,k) = std(Down_Oscillatory_amplitudeVariable(ind_down,k-1))/sqrt(length(ind_down));


    end


    figure
    if length(ind_up) >= 5
        errorbar([0:5:50],mean_vector(1,:),mean_vector(2,:),'o','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.19 .39 .68],'Color',[.19 .39 .68],'MarkerEdgeColor','k');
        hold on;
        errorbar([0:5:50],mean_vector(3,:),mean_vector(4,:),'o','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.12 0.56 1],'Color',[.12 0.56 1],'MarkerEdgeColor','k');
        hold on;
%         legend('Up. Duration','Up. Amplitude')
    end
    if length(ind_down) >= 5
        errorbar([0:5:50],mean_vector(5,:),mean_vector(6,:),'^','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.8 0.37 .8],'Color',[.8 0.37 .8],'MarkerEdgeColor','k');
        hold on;
        errorbar([0:5:50],mean_vector(7,:),mean_vector(8,:),'^','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[1 .5 1],'Color',[1 .5 1],'MarkerEdgeColor','k');
%         legend('Down. Duration', 'Down. Amplitude')
    end

%     if length(ind_down) >= 5 && length(ind_up) >= 5
%         legend('Up. Fixed Amplitude','Up. Fixed Duration', 'Down. Fixed Amplitude', 'Down. Fixed Duration')
%     end
    str = sprintf('Noise: [%g, %g)',i_test(n),j_test(n));
    title(str)
    xlabel('Area under the stimuli curve')
    ylabel('Average duration (s)')
    xlim([-2 55])
    ylim([0 200])
    
    set(gcf,'Position',[500 500 700 400])
    set(gca,'FontSize',18,'Box','off','LineWidth',1.5)
    sva = sprintf('Trial_%g.png',n);
%     saveas(gcf,sva)
    
end



%% Prueba

start = 44000;
duration_s = 50/0.01;
stimuli = zeros(1,100001);
stimuli(start:start + duration_s) = 0.1;


options = odeset('OutputFcn',@(t,y,flag) myOutputFcnStimuli(t,y,flag,w,b,I,new_memory(:,32,1),noise_t,x0,stimuli));
[t,ra] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,new_memory(:,32,1),noise_t,x0,k,r0,stimuli), [0:dt:duration], ii_cc,options);
r = ra(:,1);
plot(t,r,'-o');
% figure
% plot(t,r_inf(1:(end-1)),'-o');

% figure; plot(noise_t(40000:70000),r_inf(40000:70000))
figure; plot(new_memory(40000:65000,32,1),r_inf(40000:65000),'o','MarkerFaceColor',[.12 0.56 1]); hold on, 
plot(new_memory(44000:49000,32,1),r_inf(44000:49000),'o','MarkerFaceColor',[1 0.647 0])
xlim([-1.25 1.5])
legend('No pulse','Stimuli')
xlabel('noise'); ylabel('R infinity')
%%
s = zeros(1,130);
s(20:100) = 1;
plot(s,'LineWidth',10,'Color','k')
ylim([0 1.25])
% xticks([0 20 40 60 80 ])
% xticklabels({'60','80','100','120','140'})
xlim([0 130])
yticks([0 1])
yticklabels({'0','1'})
set(gca,'LineWidth',2,'Box','on','FontSize',50)
xlabel('Time (AU)'); ylabel('Amplitude')

%%

[t,ra0] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise_test,noise_t,x0,k,r0), [0:dt:duration], ii_cc);
[t,ra] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,new_memory(:,1,1),noise_t,x0,k,r0), [0:dt:duration], ii_cc);
[t,ra1] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,new_memory(:,2,1),noise_t,x0,k,r0), [0:dt:duration], ii_cc);
[t,ra2] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,new_memory(:,3,1),noise_t,x0,k,r0), [0:dt:duration], ii_cc);
r0 = ra0(:,1);
r = ra(:,1);
r1 = ra1(:,1);
r2 = ra2(:,1);
plot(t,r0,'.');
hold on; plot(t,r,'-o');
hold on; plot(t,r1,'-o');
hold on; plot(t,r2,'-o');
xlim([300 550])
% legend('Firing rate','noise')
%plot(t,r_inf(1:(end-1)),'-o');

%%

cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/test/15-05/UP'
load('ii_cc.mat')
load('new_memory.mat')

cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/codeTower'

start = 44000;
duration_s = 50/0.01;
stimuli = zeros(1,100001);
stimuli(start:start + duration_s) = 0.1;

stimuli5 = zeros(1,100001);
stimuli50 = zeros(1,100001);
stimuli5(start:start + duration_s) = 0.1;
stimuli50(start:start + duration_s) = 1;

[t,ra] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,new_memory(:,2,1),noise_t,x0,k,r0,stimuli), [0:dt:duration], ii_cc);
r_2 = ra(:,1);

[t,ra] = ode45(@(t,ra) rateAdapatationStimuli(t,ra,tau_r,tau_a,w,b,I,new_memory(:,32,1),noise_t,x0,k,r0,stimuli), [0:dt:duration], ii_cc);
r_32 = ra(:,1);


subplot(1,2,1)
plot(t(40000:65000),(new_memory(40000:65000,2,1) + stimuli50(40000:65000)'),'LineWidth',2); hold on;
plot(t(40000:65000),(new_memory(40000:65000,2,1) + stimuli5(40000:65000)'),'LineWidth',2); hold on;
plot(t(40000:65000),(new_memory(40000:65000,2,1)),'LineWidth',2); hold on;
plot(t(40000:65000),-0.75.*ones(1,25001),'--','LineWidth',2);
xlim([400 650])
xlabel('Time (AU)')
ylabel('Noise + Stimulus')
legend('AUC = 50','AUC = 5','No stimulus')

subplot(1,2,2)
plot(t(40000:65000),(new_memory(40000:65000,32,1) + stimuli50(40000:65000)'),'LineWidth',2); hold on;
plot(t(40000:65000),(new_memory(40000:65000,32,1) + stimuli5(40000:65000)'),'LineWidth',2); hold on;
plot(t(40000:65000),(new_memory(40000:65000,32,1)),'LineWidth',2); hold on;
plot(t(40000:65000),-0.75.*ones(1,25001),'--','LineWidth',2); hold on;
y1=get(gca,'ylim')
plot([458 458],y1,'--k','LineWidth',1.5)
xlim([400 650])
ylim([-1.5 1.5])
xlabel('Time (AU)')
ylabel('Noise + Stimulus')
legend('AuC = 50','AuC = 5','No stimulus')
set(gcf,'Position',[500 500 1200 400])





figure
subplot(1,2,1)
plot(t(40000:65000),r_2(40000:65000),'-o')
xlim([400 650])
ylabel('Firing rate')
xlabel('Time (AU)')
% suptitle('Noise transitions. High noise')

subplot(1,2,2)
plot(t(40000:65000),r_32(40000:65000),'-o'); hold on; plot([458 458],y1,'--k','LineWidth',1.5)
ylim([0 1])
xlim([400 650])
ylabel('Firing rate')
xlabel('Time (AU)')
% suptitle('Noise transitions. Low noise')



%%


i_up = [-0.8 -0.2 -0.1 0.05 0.15 0.2 0.25 0.35];
j_up = [-0.2 -0.1 0.05 0.15 0.2 0.25 0.35 0.6];

i_down = [-0.8 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2];
j_down = [-0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.6];

figure
rectangle('Position',[0.2 1.2 0.6 1],'Curvature',0.2,'FaceColor',[.12 0.56 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[0.8 1.2 0.1 1],'Curvature',0.2,'FaceColor',[.12 0.56 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[0.9 1.2 0.15 1],'Curvature',0.2,'FaceColor',[.12 0.56 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[1.05 1.2 0.1 1],'Curvature',0.2,'FaceColor',[.12 0.56 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[1.15 1.2 0.05 1],'Curvature',0.2,'FaceColor',[.12 0.56 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[1.2 1.2 0.05 1],'Curvature',0.2,'FaceColor',[.12 0.56 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[1.25 1.2 0.1 1],'Curvature',0.2,'FaceColor',[.12 0.56 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[1.35 1.2 0.25 1],'Curvature',0.2,'FaceColor',[.12 0.56 1],'LineWidth',2,'EdgeColor','k'); hold on;

rectangle('Position',[0.2 0 0.4 1],'Curvature',0.2,'FaceColor',[1 .5 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[0.6 0 0.1 1],'Curvature',0.2,'FaceColor',[1 .5 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[0.7 0 0.1 1],'Curvature',0.2,'FaceColor',[1 .5 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[0.8 0 0.1 1],'Curvature',0.2,'FaceColor',[1 .5 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[0.9 0 0.1 1],'Curvature',0.2,'FaceColor',[1 .5 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[1.0 0 0.1 1],'Curvature',0.2,'FaceColor',[1 .5 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[1.1 0 0.1 1],'Curvature',0.2,'FaceColor',[1 .5 1],'LineWidth',2,'EdgeColor','k'); hold on;
rectangle('Position',[1.2 0 0.4 1],'Curvature',0.2,'FaceColor',[1 .5 1],'LineWidth',2,'EdgeColor','k'); hold on;

xlim([0 2])
ylim([-0.25 2.25])

i_test = -0.8:0.1:0.5;
% j_test = -0.7:0.1:0.6;
% 
% for n = 1:14
% ind_up = find (avg_up(1,:) >= i_test(n) & avg_up(1,:) < j_test(n));
% ind_down = find (avg_down(1,:) >= i_test(n) & avg_down(1,:) < j_test(n));
% hist_up(n) = length(ind_up);
% hist_down(n) = length(ind_down);
% end

figure
subplot(2,1,1)
histogram(avg_up(1,:),10,'FaceColor',[.12 0.56 1],'LineWidth',2)
title('UP noise histogram')
xlabel('Noise')

subplot(2,1,2)
histogram(avg_down(1,:),10,'FaceColor',[1 .5 1],'LineWidth',2)
title('DOWN noise histogram')
xlabel('Noise')

%% Legend for fig9


plot(2,2,'o','LineWidth',2,'MarkerSize',20,'MarkerFaceColor',[.19 .39 .68],'Color',[.19 .39 .68],'MarkerEdgeColor','k'); hold on;
plot(2,2,'o','LineWidth',2,'MarkerSize',20,'MarkerFaceColor',[.12 0.56 1],'Color',[.12 0.56 1],'MarkerEdgeColor','k'); hold on;
% plot(2,2,'^','LineWidth',2,'MarkerSize',20,'MarkerFaceColor',[.8 0.37 .8],'Color',[.8 0.37 .8],'MarkerEdgeColor','k'); hold on;
% plot(2,2,'^','LineWidth',2,'MarkerSize',20,'MarkerFaceColor',[1 .5 1],'Color',[1 .5 1],'MarkerEdgeColor','k'); hold on;
legend('Up. Fixed Amplitude','Up. Fixed Duration', 'Down. Fixed Amplitude', 'Down. Fixed Duration')


%% Supplementary figure?


[t,ra] = ode45(@(t,ra) rateAdapatation(t,ra,tau_r,tau_a,w,b,I,noise,noise_t,x0,k,r0), [0:dt:duration], [0;0]);
r = ra(:,1);
plot(t(1:60000),r(1:60000)); hold on; plot(t(1:60000),noise(1:60000)); hold on
plot(t(1:60000),0.*ones(1,60000),'--','LineWidth',2);
ylim([-1 1])
legend('r','\xi')
% [r_thresh,r_cross,r_bihist,r_diptest] = BimodalThresh(r,'Schimdt');

%% Histograma

cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/hist/Oscillatory'
load('aux_down.mat')
load('aux_up.mat')
histogram(aux_down,10); hold on; histogram(aux_up,10);
legend('DOWN','UP')
% FF80FF Rosita
% 1E8EFF azulito


%% No stimuli for oscillatory

cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/200s_bistable'
load('Down_Oscillatory_amplitudeVariable.mat')
load('Down_Oscillatory_durationVariable.mat')
load('Down_Oscillatory_nostimuli.mat')
load('Up_Oscillatory_amplitudeVariable.mat')
load('Up_Oscillatory_durationVariable.mat')
load('Up_Oscillatory_nostimuli.mat')
load('oscillatory_avg_up_200.mat')
load('oscillatory_avg_down_200.mat')

plot(avg_up_200(1:100),Up_Oscillatory_nostimuli(1:100),'o','LineWidth',1,'MarkerEdgeColor',[.12 0.56 1],'MarkerSize',20)
hold on;
plot(avg_up(1,:)',nostimuli_up(:,1),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',20)
xlabel('noise')
ylabel('Duration (s)')
set(gca, 'YScale', 'log')%,'YMinorTick','off','YMinorGrid','off')
legend('Oscillatory','Bistable')
set(gcf,'Position',[500 500 700 400])
set(gca,'FontSize',22,'Box','off','LineWidth',1.5)
xlim([-0.8 0.6])
ylim([1 1000])
title('UP no stimuli.')

figure;
subplot(1,2,1)
plot(avg_down_200(1:100),Down_Oscillatory_nostimuli(1:100),'^','LineWidth',1,'MarkerEdgeColor',[1 0.5 1],'MarkerSize',20)
xlabel('noise')
ylabel('Duration (s)')
set(gca, 'YScale', 'log')%,'YMinorTick','off','YMinorGrid','off')
% legend('Oscillatory','Bistable')
set(gcf,'Position',[500 500 700 400])
set(gca,'FontSize',22,'Box','off','LineWidth',1.5)
xlim([-0.8 0.6])
ylim([1 1000])
title('Oscillatory. DOWN no stimuli.')

subplot(1,2,2)
plot(avg_down(1,:)',nostimuli_down(:,1),'^','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',20)
xlabel('noise')
ylabel('Duration (s)')
set(gca, 'YScale', 'log')%,'YMinorTick','off','YMinorGrid','off')
% legend('Oscillatory','Bistable')
set(gcf,'Position',[500 500 700 400])
set(gca,'FontSize',22,'Box','off','LineWidth',1.5)
xlim([-0.8 0.6])
ylim([1 1000])
title('Bistable. DOWN no stimuli.')


%% Figuras pendiente


cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/200s_Oscillatory'
load('Down_Oscillatory_amplitudeVariable.mat')
load('Down_Oscillatory_durationVariable.mat')
load('Down_Oscillatory_nostimuli.mat')
load('Up_Oscillatory_amplitudeVariable.mat')
load('Up_Oscillatory_durationVariable.mat')
load('Up_Oscillatory_nostimuli.mat')
load('oscillatory_avg_up_200.mat')
load('oscillatory_avg_down_200.mat')


i_test = [-0.8 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4];
j_test = [-0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.6];

% cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/Paper/Figures/Fig9/Oscillatory'

p2_final = zeros(4,10);
p1_final = zeros(4,10);
p4_final = zeros(4,10);
max_slopes = zeros(4,10);
AUC_min = 0;
AUC_max = 50;
dAUC = 5;
AUC = AUC_min:dAUC:AUC_max;
f = @(p,x) p(4) + (p(1)+p(4))./(1+exp(-p(2)*(x-p(3))));
x = linspace(AUC(1),AUC(end),100);
clc
for n = 1
    mean_vector = zeros(8,10); % Fila de datos y fila de error
    ind_up = find (avg_up_200(:,1) >= i_test(n) & avg_up_200(:,1) < j_test(n));
    ind_down = find (avg_down_200(:,1) >= i_test(n) & avg_down_200(:,1) < j_test(n));

    mean_vector(1,1) = mean(Up_Oscillatory_nostimuli(ind_up));
    mean_vector(2,1) = std(Up_Oscillatory_nostimuli(ind_up))/sqrt(length(ind_up));
    mean_vector(3,1) = mean_vector(1,1);
    mean_vector(4,1) = mean_vector(2,1);

    mean_vector(5,1) = mean(Down_Oscillatory_nostimuli(ind_down));
    mean_vector(6,1) = std(Down_Oscillatory_nostimuli(ind_down))/sqrt(length(ind_down));
    mean_vector(7,1) = mean_vector(5,1);
    mean_vector(8,1) = mean_vector(6,1);

    for k = 2:11
        mean_vector(1,k) = mean(Up_Oscillatory_durationVariable(ind_up,k-1));
        mean_vector(2,k) = std(Up_Oscillatory_durationVariable(ind_up,k-1))/sqrt(length(ind_up));

        mean_vector(3,k) = mean(Up_Oscillatory_amplitudeVariable(ind_up,k-1));
        mean_vector(4,k) = std(Up_Oscillatory_amplitudeVariable(ind_up,k-1))/sqrt(length(ind_up));

        mean_vector(5,k) = mean(Down_Oscillatory_durationVariable(ind_down,k-1));
        mean_vector(6,k) = std(Down_Oscillatory_durationVariable(ind_down,k-1))/sqrt(length(ind_down));

        mean_vector(7,k) = mean(Down_Oscillatory_amplitudeVariable(ind_down,k-1));
        mean_vector(8,k) = std(Down_Oscillatory_amplitudeVariable(ind_down,k-1))/sqrt(length(ind_down));


    end

    if n >= 4
        p1_o = max([mean_vector(1,1) mean_vector(1,end)]);
        p4_o = 0%min([mean_vector(1,1) mean_vector(1,end)]);
        slope = diff(mean_vector(1,:))/dAUC;
        [~, loc] = max(abs(slope));
        p2_o = 0%;slope(loc);
        p3_o = (AUC(loc)+AUC(loc+1))/2;
        
        p_fitted = nlinfit(AUC,mean_vector(1,:),f,[p1_o p2_o p3_o p4_o]);
        p2_final(1,n) = p_fitted(2);
        p1_final(1,n) = p_fitted(1);
        p4_final(1,n) = p_fitted(4);
        max_slopes(1,n) = max(abs(slope));

        
        
        p1_o = max([mean_vector(3,1) mean_vector(3,end)]);
        p4_o = 0;%min([mean_vector(3,1) mean_vector(3,end)]);
        slope = diff(mean_vector(3,:))/dAUC;
        [~, loc] = max(abs(slope));
        p2_o = 0;%slope(loc);
        p3_o = (AUC(loc)+AUC(loc+1))/2;
        
        p_fitted = nlinfit(AUC,mean_vector(3,:),f,[p1_o p2_o p3_o p4_o]);
        p2_final(2,n) = p_fitted(2);
        p1_final(2,n) = p_fitted(1);
        p4_final(2,n) = p_fitted(4);
        max_slopes(2,n) = max(abs(slope));
    end
    
    if n <= 6
        p1_o = max([mean_vector(5,1) mean_vector(5,end)]);
        p4_o = min([mean_vector(5,1) mean_vector(5,end)]);
        slope = diff(mean_vector(5,:))/dAUC;
        [~, loc] = max(abs(slope));
        p2_o = 0;%slope(loc);
        p3_o = (AUC(loc)+AUC(loc+1))/2;
        
        p_fitted = nlinfit(AUC,mean_vector(5,:),f,[p1_o p2_o p3_o p4_o]);
        p2_final(3,n) = p_fitted(2);
        p1_final(3,n) = p_fitted(1);
        p4_final(3,n) = p_fitted(4);
        max_slopes(3,n) = max(abs(slope));
        
        
        p1_o = max([mean_vector(7,1) mean_vector(7,end)]);
        p4_o = min([mean_vector(7,1) mean_vector(7,end)]);
        slope = diff(mean_vector(7,:))/dAUC;
        [~, loc] = max(abs(slope));
        p2_o = 0;%slope(loc);
        p3_o = (AUC(loc)+AUC(loc+1))/2;
        max_slopes(4,n) = max(abs(slope));
        
        p_fitted = nlinfit(AUC,mean_vector(7,:),f,[p1_o p2_o p3_o p4_o]);
        p2_final(4,n) = p_fitted(2);
        p1_final(4,n) = p_fitted(1);
        p4_final(4,n) = p_fitted(4);
    end
    
% f_final = p_fitted(1)./(1+exp(-p_fitted(2)*(x-p_fitted(3))));
plot(AUC,mean_vector(1,:),'o','MarkerSize',15,'LineWidth',1.5,'MarkerFaceColor',[.19 .39 .68],'Color',[.19 .39 .68],'MarkerEdgeColor','k'); hold on; 
plot(x,f(p_fitted,x),'-','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor','k')
str = sprintf('Oscillatory. Fitting. Noise: [%g, %g)',i_test(n),j_test(n));
title(str)
xlim([-2 55])
% ylim([0 200])
legend('Data','Fitting')
set(gcf,'Position',[500 500 700 400])
set(gca,'FontSize',18,'Box','off','LineWidth',1.5)
xlabel('AUC')
ylabel('Average duration (s)')
%     cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/Fits Incorrectos'
%     sva = sprintf('Bistable_4.png');
%     saveas(gcf,sva)
    
end



%% Oscillatory no noise

cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/Oscilatorio_frecuenciasVariable_noNoise'
load('osc_nonoise_DOWN_frecuenciasVariable.mat')
load('osc_nonoise_UP_frecuenciasVariable.mat')
load('AUCs_nonoise.mat')


% cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/codeTower'


locs = find(osc_nonoise_DOWN_frecuenciasVariable(1,:) > 0);
down_frecuenciasVariable = osc_nonoise_DOWN_frecuenciasVariable(:,locs);

locs = find(osc_nonoise_UP_frecuenciasVariable(1,:) > 0);
up_frecuenciasVariable = osc_nonoise_UP_frecuenciasVariable(:,locs);



mean_down = mean(down_frecuenciasVariable');
% mean_down_amplitude(2,:) = std(down_amplitudeVariable')/sqrt(length(down_amplitudeVariable'));

mean_up = mean(up_frecuenciasVariable');

plot([0 AUCs_nonoise],mean_up,'o','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.12 0.56 1],'Color',[.12 0.56 1],'MarkerEdgeColor','k'); hold on;
plot([0 AUCs_nonoise],mean_down,'^','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[1 .5 1],'Color',[1 .5 1],'MarkerEdgeColor','k');
set(gcf,'Position',[500 500 700 400])
set(gca,'FontSize',18,'Box','off','LineWidth',1.5)
ylim([0 100])
title('Oscillatory durations. No noise scenario')
xlabel('AUC')
ylabel('Average duration (s)')

% f = @(p,x) p(1)./(1+exp(-p(2)*(x-p(3))));
% 
% x = [0:5:50];
% p1_o = max(mean_up_duration);
% [M,I] = max(abs(diff(mean_up_duration)));
% p2_o = M/(x(I+1)-x(I));
% p3_o = 30%x(I+1);
% opts = statset('nlinfit'); opts.MaxIter = 1000;
% p_fitted = nlinfit(x,mean_up_duration,f,[p1_o p2_o p3_o],opts);
% 
% f_final = p_fitted(1)./(1+exp(-p_fitted(2)*(x-p_fitted(3))));
% figure;
% plot(x,mean_up_duration,'o'); hold on; plot(x,f_final,'o')

%% Figuras current injection

cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/Bistable_frecuenciasVariable/200s'
load('bistable_avg_up_200.mat'); 
load('bistable_avg_down_200.mat'); 

load('Down_Bistable_nostimuli.mat');
load('Up_Bistable_nostimuli.mat');

load('Down_frecuenciasVariable.mat');
load('Up_frecuenciasVariable.mat');

load('AUCs_down.mat');
load('AUCs_up.mat');

cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/Bistable_frecuenciasVariable/200s/figures'


i_up = [-0.8 -0.2 -0.1 0.05 0.15 0.2 0.25 0.35];
j_up = [-0.2 -0.1 0.05 0.15 0.2 0.25 0.35 0.6];

i_down = [-0.8 -0.4 -0.3 -0.2 -0.1 0 0.1 0.2];
j_down = [-0.4 -0.3 -0.2 -0.1 0 0.1 0.2 0.6];

i_test = [-0.8 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4];
j_test = [-0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.6];


for n = 2
    mean_vector = zeros(8,22); % Fila de datos y fila de error
    ind_up = find (avg_up_200 >= i_test(n) & avg_up_200 < j_test(n));
    ind_down = find (avg_down_200 >= i_test(n) & avg_down_200 < j_test(n));

    mean_vector(1,1) = mean(Up_Bistable_nostimuli(ind_up));
    test(:,1) = Up_Bistable_nostimuli(ind_up);
    mean_vector(2,1) = std(Up_Bistable_nostimuli(ind_up))/sqrt(length(ind_up));

    mean_vector(5,1) = mean(Down_Bistable_nostimuli(ind_down));
    mean_vector(6,1) = std(Down_Bistable_nostimuli(ind_down))/sqrt(length(ind_down));

    for k = 2:22
        test(:,k) = Up_frecuenciasVariables(ind_up,k-1);
        mean_vector(1,k) = mean(Up_frecuenciasVariables(ind_up,k-1));
        mean_vector(2,k) = std(Up_frecuenciasVariables(ind_up,k-1))/sqrt(length(ind_up));

 
        mean_vector(5,k) = mean(Down_frecuenciasVariables(ind_down,k-1));
        mean_vector(6,k) = std(Down_frecuenciasVariables(ind_down,k-1))/sqrt(length(ind_down));

    end


    figure
    errorbar(AUCs_up,mean_vector(1,:),mean_vector(2,:),'o','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.12 0.56 1],'Color',[.12 0.56 1],'MarkerEdgeColor','k');
    hold on;
    errorbar(AUCs_down,mean_vector(5,:),mean_vector(6,:),'^','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[1 .5 1],'Color',[1 .5 1],'MarkerEdgeColor','k');
    
    str = sprintf('Current injection. Bistable regime. Noise: [%g, %g)',i_test(n),j_test(n));
    title(str)
    xlabel('Area under the stimuli curve')
    ylabel('Average duration (s)')
    xlim([-2 55])
    ylim([0 200])
%     legend('Up','Down')
    set(gcf,'Position',[500 500 700 400])
    set(gca,'FontSize',18,'Box','off','LineWidth',1.5)
%     ylim([0 100])
    xlabel('AUC')
    ylabel('Average duration (s)')
    sva = sprintf('Final_%g.png',n);
%     saveas(gcf,sva)
    

    figure;
    plot(AUCs_up,test(:,:),'o','MarkerSize',20,'LineWidth',1.5);
    title('Trials [-0.4 -0.3)')
    xlabel('Area under the stimuli curve')
    ylabel('Average duration (s)')
    xlim([-2 55])
    ylim([0 200])
%     legend('Up','Down')
    set(gcf,'Position',[500 500 700 400])
    set(gca,'FontSize',18,'Box','off','LineWidth',1.5)
%     ylim([0 100])
    xlabel('AUC')
    ylabel('Average duration (s)')
    sva = sprintf('Final_%g.png',n);
end

%% p2 vs noise

% cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/200s_Bistable'
% load('final_values_bistable');
% load('max_slopes_bistable');

cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/Sensitivity/Data'
load('st_osc_sensitivity');
load('st_bis_sensitivity');

x_up = [-0.6 -0.35 -0.25 -0.15 -0.05 0.05 0.15 0.25 0.35 0.5];
x_down = [-0.6 -0.35 -0.25 -0.15 -0.05 0.05 0.15 0.25 0.35 0.5];

plot(x_up,abs(st_osc_sensitivity(1,:)),'o','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.12 0.56 1],'Color',[.12 0.56 1],'MarkerEdgeColor','k'); hold on;
plot(0,0.5877,'o','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.82 .82 .82],'Color',[.82 .82 .82],'MarkerEdgeColor','k'); hold on;
plot(x_up,abs(st_osc_sensitivity(2,:)),'^','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[1 .5 1],'Color',[1 .5 1],'MarkerEdgeColor','k'); hold on;
plot(0,2.8788,'^','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.82 .82 .82],'Color',[.82 .82 .82],'MarkerEdgeColor','k'); hold on;
xlabel('Noise Bins')
ylabel('Maximum Slope')
set(gcf,'Position',[500 500 700 400])
set(gca,'FontSize',18,'Box','off','LineWidth',1.5)
xlim([-0.8 0.6])
legend('UP','UP no noise','DOWN','DOWN no noise')
ylim([0 3.7])
title('Oscillatory. Spike Train input.')
sva = sprintf('st_down_sensitivity.png');
% saveas(gcf,sva)


figure; 
plot(x_up,abs(st_osc_sensitivity(1,:)),'o','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.12 0.56 1],'Color',[.12 0.56 1],'MarkerEdgeColor','k'); hold on;
plot(x_up,abs(st_bis_sensitivity(1,:)),'o','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.19 .39 .68],'Color',[.19 .39 .68],'MarkerEdgeColor','k'); hold on;
xlabel('Noise Bins')
ylabel('p2*(p1 - p4)/4')
set(gcf,'Position',[500 500 700 400])
set(gca,'FontSize',18,'Box','off','LineWidth',1.5)
xlim([-0.8 0.6])
ylim([0 3.5])
legend('Oscillatory','Bistable')
title('UP. Spike Train input.')
sva = sprintf('st_up_sensitivity.png');
% saveas(gcf,sva)


% figure;
% plot(x_up,abs(final_values_oscillatory(1,:)),'o','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.12 0.56 1],'Color',[.12 0.56 1],'MarkerEdgeColor','k'); hold on;
% % plot(x_up,abs(max_slopes(1,:)),'o','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[.12 0.56 1],'Color',[.19 .39 .68],'MarkerEdgeColor','k'); hold on;
% xlabel('Noise Bins')
% ylabel('p2*(p1 + p4)/4')
% set(gcf,'Position',[500 500 700 400])
% set(gca,'FontSize',18,'Box','off','LineWidth',1.5)
% xlim([-0.8 0.6])
% % ylim([-0.1 0.4])
% title('Oscillatory. UP')
% sva = sprintf('osc_up.png');
% saveas(gcf,sva)
% 
% 
% figure;
% plot(x_up,abs(final_values_oscillatory(3,:)),'^','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[1 .5 1],'Color',[1 .5 1],'MarkerEdgeColor','k'); hold on;
% % plot(x_up,abs(max_slopes(3,:)),'^','MarkerSize',20,'LineWidth',1.5,'MarkerFaceColor',[1 .5 1],'Color',[.19 .39 .68],'MarkerEdgeColor','k'); hold on;
% xlabel('Noise Bins')
% ylabel('p2*(p1 + p4)/4')
% set(gcf,'Position',[500 500 700 400])
% set(gca,'FontSize',18,'Box','off','LineWidth',1.5)
% xlim([-0.8 0.6])
% % ylim([-0.1 0.4])
% title('Oscillatory. DOWN')
% sva = sprintf('osc_down.png');
% saveas(gcf,sva)

%% Figuras pendiente spike train
cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/Oscilatorio_frecuenciasVariable/200s'
% load('Down_Bistable_amplitudeVariable.mat')
% load('Down_Bistable_durationVariable.mat')
load('Down_frecuenciasVariables.mat')
load('Down_Oscillatory_nostimuli.mat')
% load('Up_Bistable_amplitudeVariable.mat')
% load('Up_Bistable_durationVariable.mat')
load('Up_frecuenciasVariables.mat')
load('Up_Oscillatory_nostimuli.mat')
load('oscillatory_avg_up_200.mat')
load('oscillatory_avg_down_200.mat')

load('AUCs_down.mat')
load('AUCs_up.mat')

i_test = [-0.8 -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4];
j_test = [-0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.6];


p2_final = zeros(4,10);
p1_final = zeros(4,10);
p4_final = zeros(4,10);
max_slopes = zeros(4,10);
AUC_min = 0;
AUC_max = 50;
dAUC = 5;
AUC = AUC_min:dAUC:AUC_max;
f = @(p,x) p(4) + (p(1)+p(4))./(1+exp(-p(2)*(x-p(3))));
x = linspace(AUC(1),AUC(end),100);

st_osc_durVSauc = zeros(20,22);
s = [0 1 2 3 4 5 6 7 8 9];

for n = 1:10
    mean_vector = zeros(8,22); % Fila de datos y fila de error
    ind_up = find (avg_up_200(:,1) >= i_test(n) & avg_up_200(:,1) < j_test(n));
    ind_down = find (avg_down_200(:,1) >= i_test(n) & avg_down_200(:,1) < j_test(n));

    mean_vector(1,1) = mean(Up_Oscillatory_nostimuli(ind_up));
    mean_vector(2,1) = std(Up_Oscillatory_nostimuli(ind_up))/sqrt(length(ind_up));
    mean_vector(3,1) = mean_vector(1,1);
    mean_vector(4,1) = mean_vector(2,1);

    mean_vector(5,1) = mean(Down_Oscillatory_nostimuli(ind_down));
    mean_vector(6,1) = std(Down_Oscillatory_nostimuli(ind_down))/sqrt(length(ind_down));
    mean_vector(7,1) = mean_vector(5,1);
    mean_vector(8,1) = mean_vector(6,1);

    for k = 2:22
        mean_vector(1,k) = mean(Up_frecuenciasVariables(ind_up,k-1));
        mean_vector(2,k) = std(Up_frecuenciasVariables(ind_up,k-1))/sqrt(length(ind_up));

%         mean_vector(3,k) = mean(Up_Bistable_amplitudeVariable(ind_up,k-1));
%         mean_vector(4,k) = std(Up_Bistable_amplitudeVariable(ind_up,k-1))/sqrt(length(ind_up));

        mean_vector(5,k) = mean(Down_frecuenciasVariables(ind_down,k-1));
        mean_vector(6,k) = std(Down_frecuenciasVariables(ind_down,k-1))/sqrt(length(ind_down));

%         mean_vector(7,k) = mean(Down_Bistable_amplitudeVariable(ind_down,k-1));
%         mean_vector(8,k) = std(Down_Bistable_amplitudeVariable(ind_down,k-1))/sqrt(length(ind_down));


    end

    if n >= 4
        p1_o = max([mean_vector(1,1) mean_vector(1,end)]);
        p4_o = 0;%min([mean_vector(1,1) mean_vector(1,end)]);
        slope = diff(mean_vector(1,:))/dAUC;
        [~, loc] = max(abs(slope));
        p2_o = slope(loc);
        p3_o = (AUCs_up(loc)+AUCs_up(loc+1))/2;
        
        p_fitted_1 = nlinfit(AUCs_up,mean_vector(1,:),f,[p1_o p2_o p3_o p4_o]);
        p2_final(1,n) = p_fitted_1(2);
        p1_final(1,n) = p_fitted_1(1);
        p4_final(1,n) = p_fitted_1(4);
        max_slopes(1,n) = max(abs(slope));

        
        
%         p1_o = max([mean_vector(3,1) mean_vector(3,end)]);
%         p4_o = 0;%min([mean_vector(3,1) mean_vector(3,end)]);
%         slope = diff(mean_vector(3,:))/dAUC;
%         [~, loc] = max(abs(slope));
%         p2_o = 0;%slope(loc);
%         p3_o = (AUCs_up(loc)+AUCs_up(loc+1))/2;
%         
%         p_fitted = nlinfit(AUC,mean_vector(3,:),f,[p1_o p2_o p3_o p4_o]);
%         p2_final(2,n) = p_fitted(2);
%         p1_final(2,n) = p_fitted(1);
%         p4_final(2,n) = p_fitted(4);
%         max_slopes(2,n) = max(abs(slope));
    end
    
    if n <= 6
        p1_o = max([mean_vector(5,1) mean_vector(5,end)]);
        p4_o = min([mean_vector(5,1) mean_vector(5,end)]);
        slope = diff(mean_vector(5,:))/dAUC;
        [~, loc] = max(abs(slope));
        p2_o = 0;%slope(loc);
        p3_o = (AUCs_down(loc)+AUCs_down(loc+1))/2;
        
        p_fitted_2 = nlinfit(AUCs_down,mean_vector(5,:),f,[p1_o p2_o p3_o p4_o]);
        p2_final(3,n) = p_fitted_2(2);
        p1_final(3,n) = p_fitted_2(1);
        p4_final(3,n) = p_fitted_2(4);
        max_slopes(3,n) = max(abs(slope));
        
        
%         p1_o = max([mean_vector(7,1) mean_vector(7,end)]);
%         p4_o = min([mean_vector(7,1) mean_vector(7,end)]);
%         slope = diff(mean_vector(7,:))/dAUC;
%         [~, loc] = max(abs(slope));
%         p2_o = 0;%slope(loc);
%         p3_o = (AUCs_down(loc)+AUCs_down(loc+1))/2;
%         max_slopes(4,n) = max(abs(slope));
%         
%         p_fitted = nlinfit(AUC,mean_vector(7,:),f,[p1_o p2_o p3_o p4_o]);
%         p2_final(4,n) = p_fitted(2);
%         p1_final(4,n) = p_fitted(1);
%         p4_final(4,n) = p_fitted(4);
    end
    
% f_final = p_fitted(1)./(1+exp(-p_fitted(2)*(x-p_fitted(3))));
% plot(AUCs_up,mean_vector(1,:),'o','MarkerSize',15,'LineWidth',1.5,'MarkerFaceColor',[.19 .39 .68],'Color',[.19 .39 .68],'MarkerEdgeColor','k'); hold on; 
% plot(x,f(p_fitted_1,x),'-','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor','k')
% % plot(AUCs_up,mean_vector(5,:),'o','MarkerSize',15,'LineWidth',1.5,'MarkerFaceColor',[.19 .39 .68],'Color',[.19 .39 .68],'MarkerEdgeColor','k'); hold on; 
% % plot(x,f(p_fitted_2,x),'-','MarkerSize',15,'LineWidth',1.5,'MarkerEdgeColor','k')
% str = sprintf('Oscillatory. Fitting. Noise: [%g, %g)',i_test(n),j_test(n));
% title(str)
% xlim([-2 55])
% % ylim([0 200])
% legend('Data','Fitting')
% set(gcf,'Position',[500 500 700 400])
% set(gca,'FontSize',18,'Box','off','LineWidth',1.5)
% xlabel('AUC')
% ylabel('Average duration (s)')
%     cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/Fits Incorrectos'
%     sva = sprintf('Bistable_4.png');
%     saveas(gcf,sva)


    
st_osc_durVSauc(s(n)+n,:) = mean_vector(1,:);
st_osc_durVSauc(s(n)+n+1,:) = mean_vector(5,:);


end
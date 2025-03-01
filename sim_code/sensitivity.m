
%%% -------- Load data --------- %%%
clear all
cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/current injection/Sensitivity/Data/nonoise'
% load('st_bis_durVSauc.mat');
% load('st_bis_AUCs_up.mat');
% load('st_bis_AUCs_down.mat');
% load('st_osc_durVSauc.mat');
% load('st_osc_AUCs_up.mat');
% load('st_osc_AUCs_down.mat');
load('down_frecuencias.mat');
load('up_frecuencias.mat');
load('AUCs_nonoise.mat');

%% ---------- General sensitivity analysis ------ %%%

cd '/Users/martinesparzaiaizzo/Desktop/UPF Curso 2019-2020/2o trimestre/Practicas CBC/codeTower'
st_osc_sensitivity = zeros(2,10);

% for i = 1:10
    
    n = (1:2:20);
    y = st_osc_durVSauc(13,:);
    f = @(p,x) p(4) + (p(1)-p(4))./(1+exp(-p(2)*(x-p(3)))); 
    x = linspace(0,100);

    d = st_osc_AUCs_up;

    % Initial values and bounds

    p0 = zeros(1,4);
    p0(1) = max([y(1) y(end)]);
    p0(4) = min([y(1) y(end)]);
    slope = diff(y)/2.38;
    [~, loc] = max(abs(slope));
    p0(2) = slope(loc);
    p0(3) = (d(loc)+d(loc+1))/2;

    lb = [max(y) 0 0 min(y)];
    ub = [max(y) 0.3 50 min(y)];

    % Implementation

    error_fun = @(p_new)(p_new(4) + (p_new(1)-p_new(4))./(1+exp(-p_new(2)*(d-p_new(3)))))-y;

    p_fitted = lsqnonlin(error_fun,p0,lb,ub);
    figure;
    plot(d,y,'ko',x,f(p_fitted,x))
    legend('Data','Best fit')
    xlabel('AUC')
    xlim([0 55])
    ylabel('Duration (s)')
    (p_fitted(2)/4) * (p_fitted(1) - p_fitted(4))
%     st_osc_sensitivity(2,i) = (p_fitted(2)/4) * (p_fitted(1) - p_fitted(4));
    
% end











function [p1,p2,p3,p4,sens] = sigmoidRegression(clean,state)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

size_inp = size(clean);
f2 = -100:10:100;

% Initialize
p3 =zeros(size_inp(1),1); sens = zeros(size_inp(1),1); p4 = zeros(size_inp(1),1);
p1 = zeros(size_inp(1),1); p2 = zeros(size_inp(1),1);

for i = 1:size_inp(1)
    
    y = clean(i,:); 
    p0 = zeros(1,4);
    p0(1) = max([y(1) y(end)]);
    p0(4) = min([y(1) y(end)]);
    slope = diff(y)/5;
    [~, loc] = max(abs(slope));
    p0(2) = slope(loc);
    p0(3) = (f2(loc)+f2(loc+1))/2;

    lb = [0 -inf -100 0];
    ub = [Inf Inf 100 Inf];
    
    error_fun = @(p_new)(p_new(4) + (p_new(1)-p_new(4))./(1+exp(-p_new(2)*(f2-p_new(3)))))-y;

    p_fitted = lsqnonlin(error_fun,p0,lb,ub);
    
    if state == 'dw'
        if (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4)) < 0
            p3(i) = p_fitted(3);
            p4(i) = p_fitted(4);
            p1(i) = p_fitted(1);
            p2(i) = p_fitted(2);
            sens(i) = (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4));
        end
        if (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4)) >= 0
            p3(i) = NaN;
            p4(i) = NaN;
            p1(i) = NaN;
            p2(i) = NaN;
            sens(i) = NaN;
        end
    end
    
    if state == 'up'
        if (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4)) > 0
            p3(i) = p_fitted(3);
            p4(i) = p_fitted(4);
            p1(i) = p_fitted(1);
            p2(i) = p_fitted(2);
            sens(i) = (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4));
        end
        if (p_fitted(2) * ((p_fitted(1) + p_fitted(4))/4)) <= 0
            p3(i) = NaN;
            p4(i) = NaN;
            p1(i) = NaN;
            p2(i) = NaN;
            sens(i) = NaN;
        end
    end

%     

end
end


function [ pre,post ] = durationMeasurement( vec,start )
% durationMeasurament is a function that measures the pre and post duration
% of a signal with a given stimulus.
%   Detailed explanation goes here

    for a = (start + 1):-1:1
        if vec(a) > 0.5499
            break 
        end
    end


    for n = (start + 1):length(vec)
        if vec(n) > 0.5499
            break 
        end
    end

    pre = start + 1 - a;
    post = n - (start + 1);


end


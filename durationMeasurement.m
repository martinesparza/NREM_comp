function [ pre,post ] = durationMeasurement( vec,start, option)
% durationMeasurament is a function that measures the pre and post duration
% of a signal with a given stimulus.
%   Detailed explanation goes here

    for a = (start + 1):-1:1
        if strcmp(option,'up')
            if vec(a) < 0.5499
                break 
            end
        else
            if vec(a) > 0.5499
                break 
            end
        end
            
    end


    for n = (start + 1):length(vec)
        if strcmp(option, 'up')
            if vec(n) < 0.5499
                break 
            end
        else 
            if vec(n) > 0.5499
                break 
            end
        end

    end

    pre = start + 1 - a;
    post = n - (start + 1);


end


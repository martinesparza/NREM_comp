function [clean_out,clean_idx,fractions, idxs ] = clearCases(inp, size_inp ,stim_dur)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


over_idx = find(inp(:,ceil(size_inp(2)/2)) >= stim_dur);
remove_overstimulus = inp(over_idx,:);
over_fraction = (size_inp(1) - length(over_idx)) / size_inp(1);

sham_idx = [];
for i = 1:size_inp(1)
    if abs(inp(i,size_inp(2)) - inp(i,1)) >= 10
        sham_idx = [sham_idx i];
        
    end
end

remove_sham = inp(sham_idx,:);
sham_fraction = (size_inp(1) - length(sham_idx)) / size_inp(1);

clean_idx = [];

for i = over_idx'
    if abs(inp(i,size_inp(2)) - inp(i,1)) >= 10
        clean_idx = [clean_idx i];
        
    end
end

clean_out = inp(clean_idx,:);
clean_fraction = (size_inp(1) - length(clean_idx)) / size_inp(1);

fractions = struct('Over_Fraction',over_fraction,'Sham_Fraction',sham_fraction,...
    'Clean_Fraction', clean_fraction);

idxs = struct('Over_Idx',over_idx,'Sham_Idx',sham_idx',...
    'Clean_Idx', clean_idx');


end


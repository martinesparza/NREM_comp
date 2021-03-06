function [] = plotRaster(spikeMat, tVec,LW,color)
hold all;
for trialCount = 1:size(spikeMat,1)
    spikePos = tVec(spikeMat(trialCount, :));
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'Color',color,'LineWidth',LW);
    end
end
ylim([0 size(spikeMat, 1)+1]);
end


[spikeMat, tVec] = poissonSpikeGen(10, 10, 5,0.01);
plotRaster(spikeMat, tVec*1000);
xlabel('Time (ms)');
ylabel('Trial Number');

spikeMat = double(spikeMat);

alfa = 1/0.1;
alfa_function = alfa^2 .* tVec(1:250) .* exp(-alfa.*tVec(1:250));
figure
plot(tVec(1:250)*1000,alfa_function);

current_rate = conv(spikeMat(3,:),alfa_function);
figure;
plot(current_rate/max(current_rate))
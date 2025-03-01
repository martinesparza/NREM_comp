%% Permutation test

%% UP
mat = [sens_b; sens_o];
iterations = 10000;

test_statistic = mean(mat(1:191,1)) - mean(mat(192:end,1));
tmp_statistic = zeros(1,iterations+1);
tmp_statistic(1) = test_statistic;

for i = 2:iterations
    perm = randperm(382);
    tmp_mat = mat(perm);
    tmp_statistic(i) = mean(tmp_mat(1:191,1)) - mean(tmp_mat(192:end,1));
end

idx = find(tmp_statistic <= test_statistic);
p_value_up = length(idx)/iterations

% hist(tmp_statistic,100)

%% DOWN

mat = [p3_values_b; p3_values_o];
iterations = 10000;

test_statistic = mean(mat(1:229,1)) - mean(mat(230:end,1));
tmp_statistic = zeros(1,iterations+1);
tmp_statistic(1) = test_statistic;

for i = 2:iterations
    perm = randperm(458);
    tmp_mat = mat(perm);
    tmp_statistic(i) = mean(tmp_mat(1:229,1)) - mean(tmp_mat(230:end,1));
end

idx = find(tmp_statistic <= test_statistic);
p_value_down = length(idx)/iterations;

hist(tmp_statistic,100)
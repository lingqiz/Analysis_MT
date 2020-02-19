%% Single neuron example
nNeuron = 470;
nParas  = 5;

idx = randi(nNeuron);
parameter = fitPara(idx, :);

tuning = @(stim) tuningGauss(parameter(1), parameter(2), parameter(3), parameter(4), parameter(5), stim);

xRange = 0.1 : 0.01 : 50;
[fx, dfdx] = tuning(xRange);
numDiff = gradient(fx, xRange);

figure(); subplot(1, 2, 1); hold on;
plot(xRange, fx, '-k', 'LineWidth', 2);
plot(xRange, dfdx, '-k', 'LineWidth', 2);
plot(xRange, numDiff, '--r', 'LineWidth', 2);


transformed = log(xRange + 0.3);
numDiff = gradient(fx, transformed);

subplot(1, 2, 2); hold on;
plot(transformed, fx, '-k', 'LineWidth', 2);
plot(transformed, numDiff, '--r', 'LineWidth', 2);

%% Analysis for all the neurons
xRange = 0.1 : 0.01 : 50;
transformed = log(xRange + 1);

totalFisher = zeros(1, length(xRange));
totalTrans  = zeros(1, length(xRange));

for idx = 1:nNeuron
    parameter = fitPara(idx, :);
    tuning = @(stim) tuningGauss(parameter(1), parameter(2), parameter(3), parameter(4), parameter(5), stim);
    
    [fx, dfdx] = tuning(xRange);
    fisher = abs(dfdx) ./ sqrt(fx);
    totalFisher = totalFisher + fisher .^ 2;
    
    numDiff = gradient(fx, transformed);
    fisher = abs(numDiff) ./ sqrt(fx);
    totalTrans = totalTrans + fisher .^ 2;
end
totalFisher = sqrt(totalFisher);
totalTrans  = sqrt(totalTrans);

figure(); subplot(1, 2, 1);
plot(xRange, totalFisher, 'k', 'LineWidth', 2);

subplot(1, 2, 2);
plot(log(xRange), log(totalFisher), 'k', 'LineWidth', 2);

figure();
plot(transformed, totalTrans ./ trapz(transformed, totalTrans), 'k', 'LineWidth', 2);
ylim([0, 1]);

%% CDF transformation
transformed = cumtrapz(xRange, totalFisher);

figure();
plot(xRange, transformed, '-k', 'LineWidth', 2);

for idx = 1 : nNeuron
    parameter = fitPara(idx, :);
    tuning = @(stim) tuningGauss(parameter(1), parameter(2), parameter(3), parameter(4), parameter(5), stim);

    numDiff = gradient(fx, transformed);
    fisher = abs(numDiff) ./ sqrt(fx);
    totalTrans = totalTrans + fisher;
end

figure();
plot(transformed, totalTrans, 'k', 'LineWidth', 2);

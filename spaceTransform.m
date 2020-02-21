%% Single neuron example
nNeuron = 470;
load('./fitPara_gauss.mat');

idx = randi(nNeuron);
parameter = fitPara(idx, :);

tuning = @(stim) tuningGauss(parameter(1), parameter(2), parameter(3), parameter(4), parameter(5), stim);

xRange = 0.1 : 0.01 : 50;
[fx, dfdx] = tuning(xRange);
fisher = (dfdx .^ 2) ./ fx;

numDiff = gradient(fx, xRange);

figure(); subplot(1, 2, 1); hold on;
plot(xRange, fx, '-k', 'LineWidth', 2);
plot(xRange, dfdx, '-k', 'LineWidth', 2);
plot(xRange, numDiff, '--r', 'LineWidth', 2);

subplot(1, 2, 2);
plot(xRange, fisher, '-k', 'LineWidth', 2);

% transformed = log(xRange + 0.3);
transformed = cumtrapz(xRange, sqrt(fisher));
numDiff = gradient(fx, transformed);
fisher = (numDiff .^ 2) ./ fx;

figure();
subplot(1, 2, 1); hold on;
plot(transformed, fx, '-k', 'LineWidth', 2);
plot(transformed, numDiff, '--r', 'LineWidth', 2);

subplot(1, 2, 2);
plot(transformed, fisher, '-k', 'LineWidth', 2);

%% Analysis for all the neurons
nNeuron = 470;
load('./fitPara_gauss.mat');

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
transformed = cumtrapz(xRange, totalFisher ./ trapz(xRange, totalFisher));
numDiff = gradient(xRange, transformed);

figure();
plot(xRange, totalFisher .* numDiff, '-k', 'LineWidth', 2);
ylim([350, 450]);

totalTrans  = zeros(1, length(xRange));
for idx = 1 : nNeuron
    parameter = fitPara(idx, :);
    tuning = @(stim) tuningGauss(parameter(1), parameter(2), parameter(3), parameter(4), parameter(5), stim);
    
    fx = tuning(xRange);
    numDiff = gradient(fx, transformed);
    fisher = abs(numDiff) ./ sqrt(fx);
    totalTrans = totalTrans + fisher .^ 2;
end
totalTrans = sqrt(totalTrans);

figure();
plot(transformed, totalTrans ./ (trapz(transformed, totalTrans)), 'k', 'LineWidth', 2);
xlim([0, 1]);
ylim([0.5, 1.5]);


%% Prior transformation
load('CombinedFit/combinedMapping.mat');
prior = priorHandle(paraSub);
transformed = cumtrapz(xRange, prior(xRange));

totalTrans  = zeros(1, length(xRange));
for idx = 1 : nNeuron
    parameter = fitPara(idx, :);
    tuning = @(stim) tuningGauss(parameter(1), parameter(2), parameter(3), parameter(4), parameter(5), stim);
    
    fx = tuning(xRange);
    numDiff = gradient(fx, transformed);
    fisher = abs(numDiff) ./ sqrt(fx);
    totalTrans = totalTrans + fisher .^ 2;
end
totalTrans = sqrt(totalTrans);

figure();
plot(transformed, totalTrans ./ (trapz(transformed, totalTrans)), 'k', 'LineWidth', 2);
xlim([0, 1]);
ylim([0.0, 2.0]);

set(gca,'box','off');
set(gca,'TickDir','out');
xticks(0 : 0.2 : 1); xlabel('v tilde');
yticks(0 : 0.4 : 2.0); ylabel('fisher info');

%% Helper function
function prior = priorHandle(para)

c0 = para(1); c1 = para(2); c2 = para(3);
domain = 0.1 : 0.01 : 50;

priorUnm  = 1.0 ./ ((abs(domain) .^ c0) + c1) + c2;
nrmConst  = 1.0 / (trapz(domain, priorUnm));
prior = @(support) (1.0 ./ ((abs(support) .^ c0) + c1) + c2) * nrmConst;

end

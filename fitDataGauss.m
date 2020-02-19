%% Fit to a single neuron
neurData = load('./DeAngelis/m1c125r2.mat');

figure();
scatter(neurData.speed_values, neurData.response_values, 'k');

[parameter, func] = fitGauss(neurData.speed_values, neurData.response_values, 'rmse', 'fminsearch');

axisLim = xlim;
xRange = axisLim(1) : 0.05 : axisLim(2);

hold on;
plot(xRange, func(xRange), 'k', 'LineWidth', 1);

% Differentiation
[fx, dfdx] = func(xRange);

figure();
plot(xRange, fx, '-k'); hold on;
plot(xRange, dfdx, '--k');

% Fisher information
figure();
fisher = abs(dfdx) ./ sqrt(fx);
plot(xRange, fisher, '-k');

%% Fit to the entire population of neurons
nNeuron = 470;
nParas  = 5;
showPlot = false;

fitRSquared = zeros(nNeuron, 1);
rSquaredAvg = zeros(nNeuron, 1);
fitPara = zeros(nNeuron, nParas);

files = dir('./SpeedData/*.mat');
idx = 1;
for file = files'
    neurData = load(fullfile('./SpeedData', file.name));
    
    if showPlot
        figure();
        scatter(neurData.speed_values, neurData.response_values, 'k');
    end
    
    [parameter, func, rSquared] = fitGauss(neurData.speed_values, neurData.response_values, 'rmse', 'fminsearch');
    fitRSquared(idx) = rSquared;
    fitPara(idx, :) = parameter;
    
    speed = unique(neurData.speed_values);
    response = [];
    for vid = 1:length(speed)
        response = [response mean(neurData.response_values(neurData.speed_values == speed(vid)))];
    end
    
    sTotal = sum((response - mean(response)) .^ 2);
    sRes   = sum((response' - func(speed)) .^ 2);
    rSquaredAvg(idx) = 1 - sRes / sTotal;
    
    if showPlot
        axisLim = xlim;
        xRange = axisLim(1) : 0.05 : axisLim(2);
        hold on;
        plot(xRange, func(xRange), 'k', 'LineWidth', 1);
    end
    
    idx = idx + 1;
end

%% Goodness-of-fit/R-squared
figure(); subplot(1, 2, 1);
histogram(fitRSquared); hold on;
plot(median(fitRSquared) * ones(1, 2), ylim(), '--k', 'LineWidth', 2);
xlim([0, 1]); title(strcat('R squared:', num2str(median(fitRSquared))));

subplot(1, 2, 2);
histogram(rSquaredAvg); hold on;
plot(median(rSquaredAvg) * ones(1, 2), ylim(), '--k', 'LineWidth', 2);
xlim([0, 1]); title(strcat('R squared:', num2str(median(rSquaredAvg))));

% Scatter plot of parameter
figure();
scatter(fitPara(:, 5), fitPara(:, 3), 'MarkerEdgeColor', 'black');
xlim([0, 50]); ylim([0, 4]);

figure();
scatter(fitPara(:, 5), fitPara(:, 4), 'MarkerEdgeColor', 'black');
xlim([0, 50]); ylim([0, 3]);

%% Fisher information
showPlot = false;

xRange = 0.1 : 0.01 : 100;
totalFisher = zeros(1, length(xRange));

for idx = 1 : nNeuron
    parameter = fitPara(idx, :);
    tuning = @(stim) tuningGauss(parameter(1), parameter(2), parameter(3), parameter(4), parameter(5), stim);
    
    % Fisher information
    [fx, dfdx] = tuning(xRange);
    fisher = abs(dfdx) ./ sqrt(fx);
    
    totalFisher = totalFisher + fisher .^ 2;
    if showPlot
        figure(); subplot(1, 2, 1);
        plot(xRange, fx, '-k', 'LineWidth', 2); hold on;
        plot(xRange, dfdx, '--k', 'LineWidth', 2);
        xlim([0, 20]);
        
        subplot(1, 2, 2);
        plot(xRange, fisher, '-k', 'LineWidth', 2);
        xlim([0, 20]);
    end
end
totalFisher = sqrt(totalFisher);

normcst = trapz(xRange, totalFisher) * 2;
totalFisher = totalFisher / normcst;

% Fisher information
figure(); subplot(1, 2, 1);
plot(xRange, totalFisher, '-k', 'LineWidth', 2);

subplot(1, 2, 2);
plot(log(xRange), log(totalFisher), '-k', 'LineWidth', 2);

fitlm(log(xRange'), log(totalFisher'))
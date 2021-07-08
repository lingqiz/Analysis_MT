% Setup Plotting
plotlabOBJ = plotlab();
plotlabOBJ.applyRecipe(...
    'figureWidthInches', 20, ...
    'figureHeightInches', 8);

%% Figure 1 - Piror
load('./MappingFit/new_para_map_fit/new_para_Feb9.mat');

% figure; hold on;
% set(gca, 'FontSize', 14);
subplot(1, 2, 1); hold on;
colors = get(gca,'colororder');

nSub = 4;
allPara = [paraSub1; paraSub2; paraSub4; paraSub5];

for i = 1 : nSub
    para = allPara(i, :);
    l1 = plotPrior(para, true, false, '-', ones(1, 3) * 0.8);
end

load('CombinedFit/combinedMapping.mat');
l2 = plotPrior(paraSub, true, false, '-', ones(1, 3) * 0.1);

priorXlim = xlim();

%% Add fisher information 1
% Fisher information with Gamma tuning curve
% nNeuron = 470;
% nParas  = 5;
% load('fitPara_gamma.mat');
%
% xRange = 0.05 : 0.01 : 100;
% totalFisher = zeros(1, length(xRange));
%
% for idx = 1 : nNeuron
%     parameter = fitPara(idx, :);
%     tuning = @(stim) tuningGamma(parameter(1), parameter(2), parameter(3), parameter(4), parameter(5), stim);
%
%     % Fisher information
%     [fx, dfdx] = tuning(xRange);
%     fisher = abs(dfdx) ./ sqrt(fx);
%
%     totalFisher = totalFisher + fisher .^ 2;
% end
% totalFisher = sqrt(totalFisher);
%
% normcst = trapz(xRange, totalFisher) * 2;
% totalFisher = totalFisher / normcst;
%
% % Fisher information
% l4 = plot(log(xRange), smooth(log(totalFisher), 0.05), '-r', 'LineWidth', 2);
% fitlm(log(xRange'), log(totalFisher'))
%
% % fitlm(log(xRange'), log(totalFisher'))

%% Add fisher information 2
% Fisher information with Log-Normal tuning curve

load('fitPara_gauss.mat');
nNeuron = 470;
nParas  = 5;

xRange = 0.01 : 0.001 : 100;
totalFisher = zeros(1, length(xRange));

for idx = 1 : nNeuron
    parameter = fitPara(idx, :);
    tuning = @(stim) tuningGauss(parameter(1), parameter(2), parameter(3), parameter(4), parameter(5), stim);
    
    % Fisher information
    [fx, dfdx] = tuning(xRange);
    fisher = abs(dfdx) ./ sqrt(fx);
    
    totalFisher = totalFisher + fisher .^ 2;
end
totalFisher = sqrt(totalFisher);

normcst = trapz(xRange, totalFisher) * 2;
totalFisher = totalFisher / normcst;

% Fisher information
l5 = plot(log(xRange(xRange > 0.05 & xRange < 35)), log(totalFisher(xRange > 0.05 & xRange < 35)), 'LineWidth', 2);
fitlm(log(xRange'), log(totalFisher'))

legend([l1, l2, l5], {'Individual Subject', 'Combined', 'MT Fisher'});

% l3 = plotPrior([1, 0.33, 0], true, false, '--', ones(1, 3) * 0.1);
% legend([l1, l2, l3, l4, l5], {'Individual Subject', 'Transformation', 'Combined', 'Gamma', 'LogNormal'});

%% Carlow & Lappe Data
fovea = true;

if fovea
    subplot(1, 2, 2);
    load('Carlow_Lappe.mat');
    % data_1a(:, 2) = data_1a(:, 2) * 0.5;
    % data_1b(:, 2) = data_1b(:, 2) * 0.5;
    scatter(data_1a(:, 1), data_1a(:, 2)); hold on;
    scatter(data_1b(:, 1), data_1b(:, 2));
    
    data_1deg = [data_1a; data_1b];
    
    loss = @(para) lossFunction(data_1deg(:, 1), data_1deg(:, 2), para);
    
    options = optimoptions('fmincon','Display','iter');
    fitPara = fmincon(loss, [1, 0.1, 1e-10], [], [], [], [], [0, 0, 0], [10, 10, 1], [], options);
    
    prior = priorHandler(fitPara(1), fitPara(2), fitPara(3));
    
    xRange = 0.05 : 0.1 : 6;
    plot(xRange, prior(xRange), 'k', 'LineWidth', 2);
    
    %% Add data to the plot 1 deg
    subplot(1, 2, 1); hold on;
    scatter(log(data_1a(:, 1)), log(data_1a(:, 2)));
    scatter(log(data_1b(:, 1)), log(data_1b(:, 2)));
    
    xRange = 0.05 : 0.01 : 40;
    l6 = plot(log(xRange), log(prior(xRange)));
    
    legend([l1, l2, l5, l6], {'Individual Subject', 'Combined', 'MT Fisher', 'Carlow&Lappe'});
    
    %% 16 Deg
else
    subplot(1, 2, 2);
    load('Carlow_Lappe.mat');
    scatter(data_16(:, 1), data_16(:, 2)); hold on;
    
    loss = @(para) norm(gampdf(data_16(:, 1), para(1), para(2)) - data_16(:, 2));
    options = optimoptions('fmincon','Display','iter');
    fitPara = fmincon(loss, [1, 1], [], [], [], [], [0, 0], [1e2, 1e2], [], options);
    
    xRange = 0.05 : 0.1 : 12;
    plot(xRange, gampdf(xRange, fitPara(1), fitPara(2)), 'k', 'LineWidth', 2);
    
    %% Add data to the plot 16 deg
    subplot(1, 2, 1); hold on;
    scatter(log(data_16(:, 1)), log(data_16(:, 2))); hold on;
    
    xRange = 0.05 : 0.01 : 10;
    l6 = plot(log(xRange), log(gampdf(xRange, fitPara(1), fitPara(2))));
    
    legend([l1, l2, l5, l6], {'Individual Subject', 'Combined', 'MT Fisher', 'Carlow&Lappe'});
end

%% Helper functions
function loss = lossFunction(dataX, dataY, para)
c0 = para(1); c1 = para(2); c2 = para(3);
prior = priorHandler(c0, c1, c2);
loss  = norm(log(dataY) - log(prior(dataX)));
end

function prior = priorHandler(c0, c1, c2)
domain    = -100 : 0.01 : 100;

priorUnm  = 1.0 ./ ((abs(domain) .^ c0) + c1) + c2;
nrmConst  = 1.0 / (trapz(domain, priorUnm));
prior = @(support) (1.0 ./ ((abs(support) .^ c0) + c1) + c2) * nrmConst;
end

function line = plotPrior(para, logSpace, transform, style, lineColor)
c0 = para(1); c1 = para(2); c2 = para(3);

domain    = -100 : 0.01 : 100;

priorUnm  = 1.0 ./ ((abs(domain) .^ c0) + c1) + c2;
nrmConst  = 1.0 / (trapz(domain, priorUnm));
prior = @(support) (1.0 ./ ((abs(support) .^ c0) + c1) + c2) * nrmConst;

UB = 35; priorSupport = (0.05 : 0.001 : UB);
if logSpace
    line = plot(log(priorSupport), log(prior(priorSupport)), style, 'LineWidth', 2, 'Color', lineColor);
    
    mdl = fitlm(log(priorSupport), log(prior(priorSupport)));
    mdl.Coefficients{2, 1}
    
    labelPos = [0.05, 0.1, 0.25, 0.5, 1, 2.0, 4.0, 8.0, 20, 40];
    xticks(log(labelPos));
    xticklabels(arrayfun(@num2str, labelPos, 'UniformOutput', false));
    
    probPos = 0.01 : 0.05 : 0.3;
    yticks(log(probPos));
    yticklabels(arrayfun(@num2str, probPos, 'UniformOutput', false));
    xlim(log([0.04, 40]))
else
    priorProb = prior(priorSupport);
    if transform
        priorProb = cumtrapz(priorSupport, priorProb);
    end
    line = plot((priorSupport), (priorProb), style, 'LineWidth', 2, 'Color', lineColor);
    xlim([0.01, UB]);
end

% ylim([-7, -0.5]);
% title('Prior Across All Subjects');
xlabel('V'); ylabel('P(V)');
end

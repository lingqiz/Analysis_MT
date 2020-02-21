%% Figure 1 - Piror
load('./MappingFit/new_para_map_fit/new_para_Feb9.mat');

figure; hold on;
set(gca, 'FontSize', 14);
colors = get(gca,'colororder');

nSub = 5;
allPara = [paraSub1; paraSub2; paraSub3; paraSub4; paraSub5];

for i = 1 : nSub
    para = allPara(i, :);
    l1 = plotPrior(para, true, false, '-', ones(1, 3) * 0.8);
end

load('CombinedFit/combinedMapping.mat');
l2 = plotPrior(paraSub, true, false, '-', ones(1, 3) * 0.1);

priorXlim = xlim();
set(gcf,'Position',[0, 0, 600, 500]);
set(gca,'TickDir','out');

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
l5 = plot(log(xRange(xRange > 0.05 & xRange < 35)), log(totalFisher(xRange > 0.05 & xRange < 35)), '-b', 'LineWidth', 2);
fitlm(log(xRange'), log(totalFisher'))

legend([l1, l2, l5], {'Individual Subject', 'Combined', 'MT Fisher'});

% l3 = plotPrior([1, 0.33, 0], true, false, '--', ones(1, 3) * 0.1);
% legend([l1, l2, l3, l4, l5], {'Individual Subject', 'Transformation', 'Combined', 'Gamma', 'LogNormal'});

%% Helper functions
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
    mdl.Coefficients{2, 1};
    
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

ylim([-7, -0.5]);
title('Prior Across All Subjects');
xlabel('V'); ylabel('P(V)');
end

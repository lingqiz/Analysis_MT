%% Setup Plotting
plotlabOBJ = plotlab();
plotlabOBJ.applyRecipe(...
    'figureWidthInches', 18, ...
    'figureHeightInches', 8);

%% Plot Prior
figure; hold on;
load('./MappingFit/new_para_map_fit/new_para_Feb9.mat');

nSub = 5;
allPara = [paraSub1; paraSub2; paraSub3; paraSub4; paraSub5];

for i = [1, 2, 4, 5]
    para = allPara(i, :);
    l1 = plotPrior(para, true, false, '-', ones(1, 3) * 0.8);
end

load('CombinedFit/combinedMapping.mat');
l2 = plotPrior(paraSub, true, false, '-', ones(1, 3) * 0.1);

priorXlim = xlim();

%% Plot Threshold 1 
figure; subplot(1, 2, 1); hold on;
load('./MappingFit/new_para_map_fit/new_para_Feb9.mat');

nSub = 5;
allPara = [paraSub1; paraSub2; paraSub3; paraSub4; paraSub5];

for i = [1, 2, 4, 5]
    para = allPara(i, :);
    plotThreahold(para, true);
end

load('CombinedFit/combinedMapping.mat');
plotThreahold(paraSub, true);

labelPos = [0.25, 0.5, 1, 2.0, 4.0, 8.0, 16, 32];
xticks(log(labelPos));
xticklabels(arrayfun(@num2str, labelPos, 'UniformOutput', false));

legend({'1', '2', '3', '4', 'Com'}, 'Location', 'northeast');

grid off;
xlabel('Speed');
ylabel('Weber Fraction');

%% Plot Threshold 2
subplot(1, 2, 2); hold on;
load('./MappingFit/new_para_map_fit/new_para_Feb9.mat');

nSub = 5;
allPara = [paraSub1; paraSub2; paraSub3; paraSub4; paraSub5];

for i = [1, 2, 4, 5]
    para = allPara(i, :);
    plotThreahold(para, false);
end

load('CombinedFit/combinedMapping.mat');
plotThreahold(paraSub, false);

legend({'1', '2', '3', '4', 'Com'}, 'Location', 'northeast');

grid off;
xlabel('Speed');
ylabel('Weber Fraction');

%% Helper function
function plotThreahold(para, logSpace)
c0 = para(1); c1 = para(2); c2 = para(3);
domain    = -100 : 0.01 : 100;

priorUnm  = 1.0 ./ ((abs(domain) .^ c0) + c1) + c2;
nrmConst  = 1.0 / (trapz(domain, priorUnm));
prior = @(support) (1.0 ./ ((abs(support) .^ c0) + c1) + c2) * nrmConst;

UB = 40; priorSupport = (0.25 : 0.001 : UB);

fraction = 1 ./ prior(priorSupport) ./ priorSupport * 0.025;

if logSpace
    plot(log(priorSupport), fraction);
else
    plot(priorSupport, fraction);
end

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

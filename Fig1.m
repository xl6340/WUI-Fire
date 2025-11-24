%% Fig1B
clear; clc;

decades  = {'1990s', '2000s', '2010s', '2020s'}; nDec = numel(decades);

figure; hold on; set(gcf, 'Position', [100, 100, 600, 550]);
set(gca, 'XScale', 'log', 'YScale', 'log');
color = [38 129 182; 150 210 176; 243 215 138; 199 70 71]/255; %blue-green-yellow-red  
for i = 1:nDec
    data = readtable(sprintf('dataFig/curve/CalFire-%s.csv', decades{i}));
  
    loglog(table2array(data(:,1)), table2array(data(:,2)), 'o', 'Color', color(i,:), ...
        'MarkerFaceColor', color(i,:), 'MarkerSize', 4, 'DisplayName', decades{i});
    hold on; 
    data = readtable(sprintf('dataFig/curve/CalFire-%s-fit.csv', decades{i}));
    loglog(table2array(data(:,1)), table2array(data(:,2)), '-', 'Color', color(i,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
    hold on;
end    
ax = gca;
ax.XAxis.MinorTick = 'off'; ax.YAxis.MinorTick = 'off';
ax.FontSize = 14; 
xlabel('Fire size (km$^2$)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Probability density', 'FontSize', 16);
xticks(10 .^ (0:1:4)); xlim([10^0,10^4]);
yticks(10 .^ (-7:2:1)); ylim([10^(-7),10^1]);

data = readtable('dataFig/beta/4Decs-CalFire.csv');
beta = data.beta; betaErr = data.betaErr;
text(10^0.6, 10^-0, '$y = \alpha \cdot x^{\beta}$', 'Interpreter','latex','FontSize',20);
text(10^1.2, 10^-5, sprintf('\\beta = %.2f \\pm %.2f', beta(1), betaErr(1)), 'Color', color(i,:), ...
     'FontSize', 16, 'Interpreter', 'tex');
text(10^2.7, 10^-3.7, sprintf('\\beta = %.2f \\pm %.2f', beta(4), betaErr(4)), 'Color', color(i,:), ...
     'FontSize', 16, 'Interpreter', 'tex');


% draw bar of beta
outerAx = gca;
outerPos = outerAx.Position;
insetWidth = 0.3;
insetHeight = 0.3;
insetLeft = outerPos(1) + outerPos(3) - insetWidth;
insetBottom = outerPos(2) + outerPos(4) - insetHeight;

axInset = axes('Position', [insetLeft, insetBottom, insetWidth, insetHeight]);

b = bar(axInset, beta, 'FaceColor', 'flat', 'BarWidth', 0.6, 'BaseValue', -2);
b.CData = color;

hold(axInset, 'on');
for k = 1:numel(beta)
    errorbar(axInset, k, beta(k), betaErr(k), 'k', ...
        'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 4);
end
hold(axInset, 'off');
axInset.FontSize = 14;
axInset.XTick = 1:numel(beta);
axInset.XTickLabel = decades;
axInset.XLim = [0.5, numel(beta)+0.5];
axInset.TickDir = 'out';
axInset.Box = 'off';

axInset.YLabel.String = '$\beta$ value';
axInset.YLabel.Interpreter = 'latex';
axInset.YLabel.FontSize = 14;
axInset.YLim = [-2, -1];
axInset.YTick = -2:0.2:-1;

exportgraphics(gcf,'dataFig/Figs/Fig1B.pdf','ContentType','vector');
close(gcf);

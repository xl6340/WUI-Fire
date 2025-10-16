%% comparision among 4 datasources, all fire perimeters
clear;clc;

dataSrc = {'CalFire', 'MTBS', 'Atlas', 'FIRED'}; nSrc = numel(dataSrc);

colorL = {[167/255, 185/255, 215/255], ...  
          [250/255, 220/255, 180/255],... 
          [222/255, 163/255, 162/255], ...  
          [211/255, 211/255, 211/255]}; % lighter colors (fill)
colorD = {[87/255, 111/255, 160/255], ...
          [227/255, 184/255, 127/255], ...
          [181/255, 121/255, 121/255], ...
          [145/255, 145/255, 145/255]}; % Darker colors (edge)
Marker = {'d', '^', 'o', 's'};

N       = NaN(nSrc,1); 
beta    = NaN(nSrc,1); 
betaErr = NaN(nSrc,1); 
pBeta   = NaN(nSrc,1); 
lgC     = NaN(nSrc,1); 
plgC    = NaN(nSrc,1); 
R2      = NaN(nSrc,1); 
RMSE    = NaN(nSrc,1); 
SSE     = NaN(nSrc,1); 
logLik  = NaN(nSrc,1); 
AIC     = NaN(nSrc,1); 
BIC     = NaN(nSrc,1); 

edges = 10 .^ (-1:0.05:5); 

figure(1); hold on; set(gcf, 'Position', [100, 100, 600, 500]);
set(gca, 'XScale', 'log', 'YScale', 'log');
for i = 1:nSrc % data source
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{i}));
    size = [data.size]';
    N(i) = numel(size);
    [counts, binEdges] = histcounts(size, edges);
    binWidth = diff(binEdges);
    probDensity = counts ./ (sum(counts) * binWidth); 
    binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
  
    loglog(binCenters, probDensity, 'o', 'Color', colorD{i}, 'Marker', Marker{i}, ...
    'MarkerFaceColor', colorL{i}, 'MarkerSize', 4, 'DisplayName', dataSrc{i});
    hold on; 
    
    x = log10(binCenters(probDensity > 0)); y = log10(probDensity(probDensity > 0));
    mdl = fitlm(x, y); 

    beta(i)        = abs(mdl.Coefficients.Estimate(2));
    betaErr(i)     = mdl.Coefficients.SE(2);
    pBeta(i)       = mdl.Coefficients.pValue(2);
    lgC(i)         = mdl.Coefficients.Estimate(1);
    plgC(i)        = mdl.Coefficients.pValue(1);
    R2(i)          = mdl.Rsquared.Ordinary;
    RMSE(i)        = mdl.RMSE;
    SSE(i)         = mdl.SSE;
    logLik(i)      = mdl.LogLikelihood;
    AIC(i)         = mdl.ModelCriterion.AIC;
    BIC(i)         = mdl.ModelCriterion.BIC;

    x_fit = linspace(min(x), max(x), 100)'; y_fit = predict(mdl, x_fit);
    x_fit_orig = 10.^(x_fit); y_fit_orig = 10.^(y_fit);
    
    loglog(x_fit_orig, y_fit_orig, 'Color', colorD{i}, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    hold on;
end

ax = gca;
ax.XAxis.MinorTick = 'off'; ax.YAxis.MinorTick = 'off';
ax.FontSize = 14; 
xlabel('Fire size (km$^2$)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Probability density', 'FontSize', 16);
xticks(10 .^ (-1:1:4)); xlim([10^-1,10^4]);
yticks(10 .^ (-7:2:1)); ylim([10^(-7),10^1]);
legend('FontSize', 14, 'Box', 'off', 'Location', 'southwest');

text(gca, 0.5, 0.02, 'All fire perimeters', ...
    'Units','normalized','HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','FontSize',14,'FontWeight','bold');


% draw bar of beta
outerAx = gca;
outerPos = outerAx.Position;
insetWidth = 0.25;
insetHeight = 0.25;
insetLeft = outerPos(1) + outerPos(3) - insetWidth;
insetBottom = outerPos(2) + outerPos(4) - insetHeight;

axInset = axes('Position', [insetLeft, insetBottom, insetWidth, insetHeight]);
b = bar(axInset, abs(beta), 'FaceColor', 'flat', 'BarWidth', 0.6);
b.CData = vertcat(colorL{:});

hold(axInset, 'on');
for i = 1:4
    x = b.XEndPoints(i);
    errorbar(axInset, x, abs(beta(i)), betaErr(i), 'k', ...
        'LineWidth', 1, 'CapSize', 4);
end
hold(axInset, 'off');

axInset.XTick = 1:4;
axInset.XTickLabel = dataSrc;
axInset.XLim = [0.5, 4.5];
axInset.TickDir = 'out';
axInset.Box = 'off';
axInset.FontSize = 12;
axInset.YLabel.String = '$\beta$ value';
axInset.YLabel.Interpreter = 'latex';
axInset.YLabel.FontSize = 12;
axInset.YLim = [1, 2];
axInset.YTick = 1:0.2:2;

exportgraphics(gcf,'Figs/FigS1A.pdf','ContentType','vector');
close(gcf);

T = table( ...
    dataSrc(:), N(:), beta(:), betaErr(:), pBeta(:), ...
    lgC(:), plgC(:), R2(:), RMSE(:), SSE(:), logLik(:), AIC(:), BIC(:), ...
    'VariableNames', {'Dateset','N','Beta','Beta_SE','Beta_p', ...
                      'Log10_C','Log10_C_p','R2','RMSE','SSE','LogLik','AIC','BIC'});
writetable(T, 'Figs/Table_FigS1A.csv');
disp(T);
%% comparision among 4 datasources, 2002-2024, with same time range
clear;clc;

dataSrc = {'CalFire', 'MTBS', 'Atlas', 'FIRED'}; nSrc = numel(dataSrc);

colorL = {[167/255, 185/255, 215/255], ...  
          [250/255, 220/255, 180/255],... 
          [222/255, 163/255, 162/255], ...  
          [211/255, 211/255, 211/255]}; % lighter colors (fill)
colorD = {[87/255, 111/255, 160/255], ...
          [227/255, 184/255, 127/255], ...
          [181/255, 121/255, 121/255], ...
          [145/255, 145/255, 145/255]}; % Darker colors (edge)
Marker = {'d', '^', 'o', 's'};

N       = NaN(nSrc,1); 
beta    = NaN(nSrc,1); 
betaErr = NaN(nSrc,1); 
pBeta   = NaN(nSrc,1); 
lgC     = NaN(nSrc,1); 
plgC    = NaN(nSrc,1); 
R2      = NaN(nSrc,1); 
RMSE    = NaN(nSrc,1); 
SSE     = NaN(nSrc,1); 
logLik  = NaN(nSrc,1); 
AIC     = NaN(nSrc,1); 
BIC     = NaN(nSrc,1); 

edges = 10 .^ (-1:0.05:5); 

figure(1); hold on; set(gcf, 'Position', [100, 100, 600, 500]);
set(gca, 'XScale', 'log', 'YScale', 'log');
for i = 1:4 % data source
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{i}));
    size = [data.size]';
    year = [data.year]';
    size = size(year > 2001 & year < 2025);
    N(i) = numel(size);

    [counts, binEdges] = histcounts(size, edges);
    binWidth = diff(binEdges);
    probDensity = counts ./ (sum(counts) * binWidth); 
    binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
  
    loglog(binCenters, probDensity, 'o', 'Color', colorD{i}, 'Marker', Marker{i}, ...
    'MarkerFaceColor', colorL{i}, 'MarkerSize', 4, 'DisplayName', dataSrc{i});
    hold on; 
    
    x = log10(binCenters(probDensity > 0)); y = log10(probDensity(probDensity > 0));
    mdl = fitlm(x, y); 

    beta(i)        = abs(mdl.Coefficients.Estimate(2));
    betaErr(i)     = mdl.Coefficients.SE(2);
    pBeta(i)       = mdl.Coefficients.pValue(2);
    lgC(i)         = mdl.Coefficients.Estimate(1);
    plgC(i)        = mdl.Coefficients.pValue(1);
    R2(i)          = mdl.Rsquared.Ordinary;
    RMSE(i)        = mdl.RMSE;
    SSE(i)         = mdl.SSE;
    logLik(i)      = mdl.LogLikelihood;
    AIC(i)         = mdl.ModelCriterion.AIC;
    BIC(i)         = mdl.ModelCriterion.BIC;

    x_fit = linspace(min(x), max(x), 100)'; y_fit = predict(mdl, x_fit);
    x_fit_orig = 10.^(x_fit); y_fit_orig = 10.^(y_fit);
    
    loglog(x_fit_orig, y_fit_orig, 'Color', colorD{i}, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    hold on;
end

ax = gca;
ax.XAxis.MinorTick = 'off'; ax.YAxis.MinorTick = 'off';
ax.FontSize = 14; 
xlabel('Fire size (km$^2$)', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('Probability density', 'FontSize', 16);
xticks(10 .^ (-1:1:4)); xlim([10^-1,10^4]);
yticks(10 .^ (-7:2:1)); ylim([10^(-7),10^1]);

legend('FontSize', 14, 'Box', 'off', 'Location', 'southwest');
text(gca, 0.5, 0.02, 'Fires between 2002-2024', ...
    'Units','normalized','HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','FontSize',14,'FontWeight','bold');

% draw bar of beta
outerAx = gca;
outerPos = outerAx.Position;
insetWidth = 0.25;
insetHeight = 0.25;
insetLeft = outerPos(1) + outerPos(3) - insetWidth;
insetBottom = outerPos(2) + outerPos(4) - insetHeight;

axInset = axes('Position', [insetLeft, insetBottom, insetWidth, insetHeight]);
b = bar(axInset, abs(beta), 'FaceColor', 'flat', 'BarWidth', 0.6);
b.CData = vertcat(colorL{:});

hold(axInset, 'on');
for i = 1:4
    x = b.XEndPoints(i);
    errorbar(axInset, x, abs(beta(i)), betaErr(i), 'k', ...
        'LineWidth', 1, 'CapSize', 4);
end
hold(axInset, 'off');

axInset.XTick = 1:4;
axInset.XTickLabel = dataSrc;
axInset.XLim = [0.5, 4.5];
axInset.TickDir = 'out';
axInset.Box = 'off';
axInset.FontSize = 12;
axInset.YLabel.String = '$\beta$ value';
axInset.YLabel.Interpreter = 'latex';
axInset.YLabel.FontSize = 12;
axInset.YLim = [1, 2];
axInset.YTick = 1:0.2:2;

exportgraphics(gcf,'Figs/FigS1B.pdf','ContentType','vector');
close(gcf);

T = table( ...
    dataSrc(:), N(:), beta(:), betaErr(:), pBeta(:), ...
    lgC(:), plgC(:), R2(:), RMSE(:), SSE(:), logLik(:), AIC(:), BIC(:), ...
    'VariableNames', {'Dateset','N','Beta','Beta_SE','Beta_p', ...
                      'Log10_C','Log10_C_p','R2','RMSE','SSE','LogLik','AIC','BIC'});
writetable(T, 'Figs/Table_FigS1B.csv');
disp(T);
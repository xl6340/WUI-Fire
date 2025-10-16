%% decade beta trend for CalFire
clear;clc;

decades = {'1990s', '2000s', '2010s', '2020s'};
nDec     = numel(decades);

N       = NaN(1,nDec);
beta    = NaN(1,nDec);
betaErr = NaN(1,nDec);
pBeta   = NaN(1,nDec);
lgC     = NaN(1,nDec);
plgC    = NaN(1,nDec);
R2      = NaN(1,nDec);
RMSE    = NaN(1,nDec);
SSE     = NaN(1,nDec);
logLik  = NaN(1,nDec);
AIC     = NaN(1,nDec);
BIC     = NaN(1,nDec);

freq1000 = NaN(1, nDec);

edges = 10 .^ (-1:0.05:5); 
data = shaperead('dataPrc/firePrmt/CalFire.shp');
decade = {data.decade}'; 
sizeAll = [data.size]';

  
figure; hold on; set(gcf, 'Position', [100, 100, 600, 550]);
set(gca, 'XScale', 'log', 'YScale', 'log');
color = {[38/255 129/255 182/255],[150/255 210/255 176/255],...
         [243/255 215/255 138/255], [199/255 70/255 71/255]}; %blue-green-yellow-red  
for i = 1:nDec        
    size = sizeAll(strcmp(decade, decades{i}));
    N(1,i) = numel(size);

    if numel(size) < 10
        continue;
    end

    [counts, binEdges] = histcounts(size, edges);
    binWidth = diff(binEdges);
    probDensity = counts ./ (sum(counts) * binWidth); 
    binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
  
    loglog(binCenters, probDensity, 'o', 'Color', color{i}, ...
        'MarkerFaceColor', color{i}, 'MarkerSize', 4, 'DisplayName', decades{i});
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

    % Calculate probability density for a 1000 km^2 fire ---
    targetSize = 1000;
    log10_prob_target = predict(mdl, log10(targetSize));
    prob1000 = 10^log10_prob_target;
    binIndex = discretize(targetSize, edges); % <-- Find the bin index
    if ~isnan(binIndex)
        targetBinWidth = edges(binIndex + 1) - edges(binIndex); 
        freq1000(i) = prob1000 * targetBinWidth; 
    end
    
    x_fit = linspace(min(x), max(x), 100)'; y_fit = predict(mdl, x_fit);
    x_fit_orig = 10.^(x_fit); y_fit_orig = 10.^(y_fit);
    
    loglog(x_fit_orig, y_fit_orig, '-', 'Color',color{i}, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    hold on;
end    
ax = gca;
ax.XAxis.MinorTick = 'off'; ax.YAxis.MinorTick = 'off';
ax.FontSize = 14; 
xlabel('Fire size (km$^2$)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Probability density', 'FontSize', 16);
xticks(10 .^ (0:1:4)); xlim([10^0,10^4]);
yticks(10 .^ (-7:2:1)); ylim([10^(-7),10^1]);


% draw bar of beta
outerAx = gca;
outerPos = outerAx.Position;
insetWidth = 0.3;
insetHeight = 0.3;
insetLeft = outerPos(1) + outerPos(3) - insetWidth;
insetBottom = outerPos(2) + outerPos(4) - insetHeight;

axInset = axes('Position', [insetLeft, insetBottom, insetWidth, insetHeight]);

validMask = ~isnan(beta(:));
vals      = beta(validMask);
errs      = betaErr(validMask);
cols      = vertcat(color{validMask});
decLabels = decades(validMask);

b = bar(axInset, vals, 'FaceColor', 'flat', 'BarWidth', 0.6);
b.CData = cols;

hold(axInset, 'on');
for k = 1:numel(vals)
    errorbar(axInset, k, vals(k), errs(k), 'k', ...
        'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 4);
end
hold(axInset, 'off');

axInset.FontSize = 14;
axInset.XTick = 1:numel(vals);
axInset.XTickLabel = decLabels;
axInset.XLim = [0.5, numel(vals)+0.5];
axInset.TickDir = 'out';
axInset.Box = 'off';

axInset.YLabel.String = '$\beta$ value';
axInset.YLabel.Interpreter = 'latex';
axInset.YLabel.FontSize = 14;
axInset.YLim = [1, 2];
axInset.YTick = 1:0.2:2;

if ~exist('Figs','dir'), mkdir Figs; end
exportgraphics(gcf,'Figs/Fig1B.pdf','ContentType','vector');
close(gcf);


% Assemble and save table
T = table( ...
    decades(:), N(:), beta(:), betaErr(:), pBeta(:), ...
    lgC(:), plgC(:), R2(:), RMSE(:), SSE(:), logLik(:), AIC(:), BIC(:), ...
    'VariableNames', {'Decade','N','Beta','Beta_SE','Beta_p', ...
                      'Log10_C','Log10_C_p','R2','RMSE','SSE','LogLik','AIC','BIC'});

writetable(T, 'Figs/TableS2_Fig1B.csv');
disp(T);



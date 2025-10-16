%% sizeRange for CalFire
clear;clc;

sizeRange = {'≤1 km^2', '≥1 km^2'};  nSize = numel(sizeRange);

N       = NaN(1,nSize);
beta    = NaN(1,nSize);
betaErr = NaN(1,nSize);
pBeta   = NaN(1,nSize);
lgC     = NaN(1,nSize);
plgC    = NaN(1,nSize);
R2      = NaN(1,nSize);
RMSE    = NaN(1,nSize);
SSE     = NaN(1,nSize);
logLik  = NaN(1,nSize);
AIC     = NaN(1,nSize);
BIC     = NaN(1,nSize);

edge{1} = 10 .^ (-4:0.05:-2);
edge{2} = 10 .^ (-2:0.05:5);

data = shaperead('dataPrc/firePrmt/CalFire_fullSize.shp');
sizeAll = [data.GIS_ACRES]' .* 0.0040468564224;

figure; hold on; set(gcf, 'Position', [100, 100, 600, 550]);
set(gca, 'XScale', 'log', 'YScale', 'log');
color = {[38/255 129/255 182/255],[150/255 210/255 176/255]}; %blue-green
for i = 1:nSize
    edges = edge{i};
    if i == 1
        size = sizeAll(sizeAll > 0 & sizeAll < 0.01);
    else
        size = sizeAll(sizeAll >= 0.01);
    end

    N(1,i) = numel(size);

    if numel(size) < 10, continue; end

    [counts, binEdges] = histcounts(size, edges);
    binWidth = diff(binEdges);
    probDensity = counts ./ (sum(counts) * binWidth); 
    binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
  
    loglog(binCenters, probDensity, 'o', 'Color', color{i}, ...
        'MarkerFaceColor', color{i}, 'MarkerSize', 4, 'DisplayName', sizeRange{i});
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
    
    loglog(x_fit_orig, y_fit_orig, '-', 'Color',color{i}, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    hold on;

    xpos = 10^(mean(x_fit));  
    ypos = 10^(mean(y_fit));  
    label = sprintf('\\beta = %.2f \\pm %.2f', abs(beta(i)), betaErr(i));
    text(xpos, ypos, label, 'Color', 'k', 'FontSize', 12, 'Interpreter', 'tex');
end    
ax = gca;
ax.XAxis.MinorTick = 'off'; ax.YAxis.MinorTick = 'off';
ax.FontSize = 14; 
xlabel('Fire size (km$^2$)', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('Probability density', 'FontSize', 16);
xticks(10 .^ (-5:1:4)); xlim([10^-5,10^4]);
yticks(10 .^ (-7:2:4)); ylim([10^(-7),10^4]);

exportgraphics(gcf,'Figs/FigS2.pdf','ContentType','vector');
close(gcf);

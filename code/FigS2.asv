%% decade beta trend for 4 datasets
clear;clc;

decades = {'1990s', '2000s', '2010s', '2020s'};
dataSrc = {'CalFire', 'MTBS', 'Atlas', 'FIRED'};

nSrc     = numel(dataSrc);
nDec     = numel(decades);

N       = NaN(nSrc,nDec);
beta    = NaN(nSrc,nDec);
betaErr = NaN(nSrc,nDec);
pBeta   = NaN(nSrc,nDec);
lgC     = NaN(nSrc,nDec);
plgC    = NaN(nSrc,nDec);
R2      = NaN(nSrc,nDec);
RMSE    = NaN(nSrc,nDec);
SSE     = NaN(nSrc,nDec);
logLik  = NaN(nSrc,nDec);
AIC     = NaN(nSrc,nDec);
BIC     = NaN(nSrc,nDec);

edges = 10 .^ (-1:0.05:5); 

figure('Position',[100 100 1200 900]);
color = {[38/255 129/255 182/255],[150/255 210/255 176/255],...
         [243/255 215/255 138/255], [199/255 70/255 71/255]}; %blue-green-yellow-red  

for s = 1:nSrc
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{s}));
    decade = {data.decade}'; 
    sizeAll = [data.size]';

    subplot(2,2,s); hold on;
    set(gca, 'XScale', 'log', 'YScale', 'log');

    for i = 1:nDec        
        size = sizeAll(strcmp(decade, decades{i}));
        N(s,i) = numel(size);

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

        beta(s,i)        = abs(mdl.Coefficients.Estimate(2));
        betaErr(s,i)     = mdl.Coefficients.SE(2);
        pBeta(s,i)       = mdl.Coefficients.pValue(2);
        lgC(s,i)         = mdl.Coefficients.Estimate(1);
        plgC(s,i)        = mdl.Coefficients.pValue(1);
        R2(s,i)          = mdl.Rsquared.Ordinary;
        RMSE(s,i)        = mdl.RMSE;
        SSE(s,i)         = mdl.SSE;
        logLik(s,i)      = mdl.LogLikelihood;
        AIC(s,i)         = mdl.ModelCriterion.AIC;
        BIC(s,i)         = mdl.ModelCriterion.BIC;

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
    xticks(10 .^ (-1:1:4)); xlim([10^-1,10^4]);
    yticks(10 .^ (-7:2:1)); ylim([10^(-7),10^1]);
    legend('FontSize', 14, 'Box', 'on', 'Location', 'southwest');

    text(0.02, 0.98, dataSrc{s},'Units', 'normalized','FontSize', 16, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    
    % inset bar plot
    outerAx = gca;
    outerPos = outerAx.Position;
    insetWidth = 0.25*outerPos(3);
    insetHeight = 0.25*outerPos(4);
    insetLeft = outerPos(1) + outerPos(3) - insetWidth - 0.02;
    insetBottom = outerPos(2) + outerPos(4) - insetHeight - 0.02;
    
    axInset = axes('Position', [insetLeft, insetBottom, insetWidth, insetHeight]);

    validMask = ~isnan(beta(s,:));
    vals = beta(s,validMask);
    errs = betaErr(s,validMask);
    cols = vertcat(color{validMask});
    decLabels = decades(validMask);

    b = bar(axInset, vals, 'FaceColor', 'flat', 'BarWidth', 0.6); b.EdgeColor = 'none';
    b.CData = cols;
    
    hold(axInset, 'on');
    for k = 1:numel(vals)
        errorbar(axInset, k, vals(k), errs(k), 'k','LineStyle', 'none', 'LineWidth', 1, 'CapSize', 4);
    end
    hold(axInset, 'off');
    
    axInset.XTick = 1:numel(vals);
    axInset.XTickLabel = decLabels;
    axInset.XLim = [0.5, numel(vals)+0.5];
    axInset.TickDir = 'out';
    axInset.Box = 'off';
    axInset.FontSize = 12;
    axInset.YLabel.String = '$\beta$ value';
    axInset.YLabel.Interpreter = 'latex';
    axInset.YLabel.FontSize = 14;
    axInset.YLim = [1, 2];
    axInset.YTick = 1:0.5:2;
end

exportgraphics(gcf,'Figs/FigS2.pdf','ContentType','vector');
close(gcf);

% export data
decades_vec = repelem(decades(:), nSrc, 1);
dataSrc_vec = repmat(dataSrc(:), nDec, 1);

T = table( ...
    decades_vec, dataSrc_vec, ...
    N(:), beta(:), betaErr(:), pBeta(:), lgC(:), plgC(:), R2(:), RMSE(:), SSE(:), logLik(:), AIC(:), BIC(:),...
    'VariableNames', {'decades','dataSrc','N','Beta','Beta_SE','Beta_p', ...
                      'Log10_C','Log10_C_p','R2','RMSE','SSE','LogLik','AIC','BIC'});
writetable(T, 'Figs/TableX_FigS2.csv');
disp(T);
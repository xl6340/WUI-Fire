clear;clc;

dataSrc = {'CalFire', 'MTBS', 'Atlas', 'FIRED'}; nSrc = numel(dataSrc);
fires = {'WUI', 'Wildland'};           nFire = numel(fires);

color = {[216, 118, 89]/255; [41, 157, 143]/255}; % WUI, Wildland

N       = NaN(nSrc, nFire);
beta    = NaN(nSrc, nFire);
betaErr = NaN(nSrc, nFire);
pBeta   = NaN(nSrc, nFire);
lgC     = NaN(nSrc, nFire);
plgC    = NaN(nSrc, nFire);
R2      = NaN(nSrc, nFire);
RMSE    = NaN(nSrc, nFire);
SSE     = NaN(nSrc, nFire);
logLik  = NaN(nSrc, nFire);
AIC     = NaN(nSrc, nFire);
BIC     = NaN(nSrc, nFire);

edges = 10 .^ (-1:0.05:5);

figure; clf; set(gcf, 'Position', [100, 100, 400, 300]); %[100, 100, 600, 450]
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); 
for s = 1:nSrc
    ax(s) = nexttile; hold on;

    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{s}));
    sizeAll = [data.size]';
    WUIFire = [data.BrdWUIFire]'; 

    for j = 1:nFire
        size = sizeAll(WUIFire == 2-j);
        N(s,j) = numel(size);

        if numel(size) < 10, continue; end

        [counts, binEdges] = histcounts(size, edges);
        binWidth = diff(binEdges);
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));

        loglog(binCenters, probDensity, 'o', 'Color', color{j},...
               'MarkerFaceColor', color{j}, 'MarkerSize', 2, 'DisplayName', fires{j},'LineStyle', 'none'); 
        
        x = log10(binCenters(probDensity > 0)); 
        y = log10(probDensity(probDensity > 0));
        mdl = fitlm(x, y); 

        beta(s,j)     = abs(mdl.Coefficients.Estimate(2));
        betaErr(s,j)  = mdl.Coefficients.SE(2);
        pBeta(s,j)    = mdl.Coefficients.pValue(2);
        lgC(s,j)      = mdl.Coefficients.Estimate(1);
        plgC(s,j)     = mdl.Coefficients.pValue(1);
        R2(s,j)       = mdl.Rsquared.Ordinary;
        RMSE(s,j)     = mdl.RMSE;
        SSE(s,j)      = mdl.SSE;
        logLik(s,j)   = mdl.LogLikelihood;
        AIC(s,j)      = mdl.ModelCriterion.AIC;
        BIC(s,j)      = mdl.ModelCriterion.BIC;

        x_fit = linspace(min(x), max(x), 100)';
        y_fit = predict(mdl, x_fit);
        loglog(10.^x_fit, 10.^y_fit, '-', 'Color', color{j}, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end

    set(gca, 'XScale', 'log', 'YScale', 'log', ...
             'FontSize', 12, 'XMinorTick', 'off', 'YMinorTick', 'off', ...
             'TickLength', [0.02 0.02], 'LineWidth', 0.25); 

    xticks(10 .^ (-1:2:4)); xlim([10^-1, 10^4]);
    yticks(10 .^ (-7:3:1)); ylim([10^(-7), 10^1]);    
    text(0.1, 0.1, dataSrc{s}, 'Units', 'normalized','HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom','FontSize', 14, 'FontWeight', 'normal');
end

set(ax(1), 'XTickLabel', []);
set(ax(2), 'XTickLabel', []);

set(ax(2), 'YTickLabel', []);
set(ax(4), 'YTickLabel', []);

xlabel(t, 'Fire size (km^2)', 'FontSize', 14);
ylabel(t, 'Probability density', 'FontSize', 14);

% exportgraphics(gcf,'Figs/Fig2A.pdf','ContentType','vector');
% close(gcf);
%% bar of beta
figure; clf; hold on;
set(gcf, 'Position', [100, 100, 400, 300]);

b = bar(beta, 'grouped', 'BarWidth', 0.8); 
for j = 1:nFire
    b(j).FaceColor = 'flat';
    b(j).EdgeColor   = 'none';
    b(j).DisplayName = fires{j};

    cd = nan(nSrc,3);
    for i = 1:nSrc, cd(i,:) = color{j}; end
    b(j).CData = cd;

    drawnow;                              % ensure XEndPoints exists
    xj = b(j).XEndPoints(:);              % nSrc×1
    yj = beta(:,j);                       % nSrc×1
    ej = betaErr(:,j);                    % nSrc×1

    mask = ~(isnan(xj) | isnan(yj) | isnan(ej));  % nSrc×1 logical
    errorbar(xj(mask), yj(mask), ej(mask), 'k', ...
        'linestyle','none','linewidth',1,'capsize',4,'HandleVisibility','off');
end

set(gca, 'FontSize', 14, 'LineWidth', 0.25);
xticks(1:nSrc);
xticklabels({'CalFire','MTBS','Atlas','FIRED'});
set(gca, 'XTickLabelRotation', 30);

ylim(1:2);
yticks(1:0.2:2);
ylabel('$\beta$ value', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

lgd = legend('Location', 'northwest', 'Box', 'off', 'FontSize', 14);
pos = lgd.Position;   % [x y width height]
pos(2) = pos(2) + 0.05;  % shift upward (adjust 0.05 as needed)
lgd.Position = pos;

% exportgraphics(gcf,'Figs/Fig2B.pdf','ContentType','vector');
% close(gcf);

% export data
Fire = {'WUI', 'WUI','WUI','WUI','Wildland','Wildland','Wildland','Wildland'};  
Src = [dataSrc dataSrc];
T = table( ...
    Fire(:), Src(:), N(:), beta(:), betaErr(:), pBeta(:), ...
    lgC(:), plgC(:), R2(:), RMSE(:), SSE(:), logLik(:), AIC(:), BIC(:), ...
    'VariableNames', {'Fire','dataSrc','N','Beta','Beta_SE','Beta_p', ...
                      'Log10_C','Log10_C_p','R2','RMSE','SSE','LogLik','AIC','BIC'});

writetable(T, 'Figs/TableS3_Fig2.csv');
disp(T);


%% decade beta trend for WUI and nonWUI
clear;clc;

dataSrc = {'CalFire', 'MTBS', 'Atlas', 'FIRED'}; nSrc = numel(dataSrc);
decades = {'1990s', '2000s', '2010s', '2020s'};  nDec = numel(decades);
fires   = {'WUI', 'Wildland'};                   nFire = numel(fires);

N       = NaN(nSrc, nDec, nFire);
beta    = NaN(nSrc, nDec, nFire);
betaErr = NaN(nSrc, nDec, nFire);
pBeta   = NaN(nSrc, nDec, nFire);
lgC     = NaN(nSrc, nDec, nFire);
plgC    = NaN(nSrc, nDec, nFire);
R2      = NaN(nSrc, nDec, nFire);
RMSE    = NaN(nSrc, nDec, nFire);
SSE     = NaN(nSrc, nDec, nFire);
logLik  = NaN(nSrc, nDec, nFire);
AIC     = NaN(nSrc, nDec, nFire);
BIC     = NaN(nSrc, nDec, nFire);

edges = 10 .^ (-1:0.05:5);
for s = 1:nSrc
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{s}));

    sizeAll = [data.size]';
    Dec = {data.decade}'; 
    WUIFire = [data.BrdWUIFire]'; 

    for i = 1:nDec
        for j = 1:nFire
            size = sizeAll(strcmp(Dec, decades{i}) & WUIFire == 2-j);
            N(s,i,j) = numel(size);
    
            if numel(size) < 10, continue; end

            [counts, binEdges] = histcounts(size, edges);
            binWidth = diff(binEdges);
            probDensity = counts ./ (sum(counts) * binWidth); 
            binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
            
            x = log10(binCenters(probDensity > 0)); 
            y = log10(probDensity(probDensity > 0));
            mdl = fitlm(x, y); 

            beta(s,i,j)     = abs(mdl.Coefficients.Estimate(2));
            betaErr(s,i,j)  = mdl.Coefficients.SE(2);
            pBeta(s,i,j)    = mdl.Coefficients.pValue(2);
            lgC(s,i,j)      = mdl.Coefficients.Estimate(1);
            plgC(s,i,j)     = mdl.Coefficients.pValue(1);
            R2(s,i,j)       = mdl.Rsquared.Ordinary;
            RMSE(s,i,j)     = mdl.RMSE;
            SSE(s,i,j)      = mdl.SSE;
            logLik(s,i,j)   = mdl.LogLikelihood;
            AIC(s,i,j)      = mdl.ModelCriterion.AIC;
            BIC(s,i,j)      = mdl.ModelCriterion.BIC;
        end
    end  
end

% Plotting trend line
color = {[0.769, 0.635, 0.820],[0.976, 0.686, 0.545], [0.314, 0.459, 0.314],[0.396, 0.690, 0.910]}; 
Marker = {'d', 'o', '^', 's'};

figure; hold on; set(gcf, 'Position', [100, 100, 550, 400]);
for s = 1:nSrc
    for j = 1:nFire
        y = squeeze(beta(s,:, j));
        e = squeeze(betaErr(s,:, j));
        
        validIdx = ~isnan(y);  
        x = 1:nDec; 
        x = x(validIdx);
        y = y(validIdx);
        e = e(validIdx);

        errorbar(x, y, e, 'Color', color{s}, 'LineStyle', 'none', ...
            'LineWidth', 1.5, 'CapSize', 4, 'HandleVisibility', 'off');

        if j ==1
            plot(x, y, '-o', 'Color', color{s}, 'MarkerFaceColor', color{s}, 'MarkerEdgeColor', 'k', ...
                'Marker', Marker{s}, 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', fires{j});
        elseif j == 2
            plot(x, y, '--o', 'Color', color{s}, 'MarkerFaceColor', color{s}, 'MarkerEdgeColor','k', ...
                'Marker', Marker{s}, 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', fires{j});
        end
    end
end
text(gca, 0.02, 0.05, 'WholeCal','Units','normalized', 'HorizontalAlignment','left', ... 
    'VerticalAlignment','bottom', 'FontSize',12, 'FontWeight','bold', ... 
    'BackgroundColor','w', 'Margin',2, 'Clipping','on');
set(gca, 'FontSize', 12);
xlim([0.8, 4.2]); 
ylim([0.95, 2]); 
% ylim([0, 2]); 
xticks(1:4); 
xticklabels(decades);
yticks(1:0.2:2);
ylabel('\beta value', 'Interpreter', 'tex', 'FontSize', 14);
set(gca, 'FontSize', 16);
box on;
legend off;

exportgraphics(gcf, 'Figs/Fig3A.pdf', 'ContentType', 'vector');
close(gcf);

% export data: NaN(nSrc, nDec, nFire);
dataSrc_vec = repmat(dataSrc(:), nDec*nFire, 1);                 % 4×(4*2) -> 32×1
decade_vec  = repmat(repelem(decades(:), nSrc, 1), nFire, 1);    % (each decade repeated nSrc) repeated for nFire
fire_vec    = repelem(fires(:), nSrc*nDec, 1);                   % each fire repeated nSrc*nDec
T = table( ...
    dataSrc_vec, decade_vec, fire_vec, ...
    N(:), beta(:), betaErr(:), pBeta(:), lgC(:), plgC(:), R2(:), RMSE(:), SSE(:), logLik(:), AIC(:), BIC(:), ...
    'VariableNames', {'dataSrc','decade','fire','N','Beta','Beta_SE','Beta_p', ...
                      'Log10_C','Log10_C_p','R2','RMSE','SSE','LogLik','AIC','BIC'});
writetable(T, 'Figs/TableS4_Fig3A.csv');
disp(T);
%% NorCal
clear;clc;

dataSrc = {'CalFire', 'MTBS', 'Atlas', 'FIRED'}; nSrc = numel(dataSrc);
decades = {'1990s', '2000s', '2010s', '2020s'}; nDec = numel(decades);
fires = {'WUI', 'Wildland'};                    nFire = numel(fires);

N       = NaN(nSrc, nDec, nFire);
beta    = NaN(nSrc, nDec, nFire);
betaErr = NaN(nSrc, nDec, nFire);
pBeta   = NaN(nSrc, nDec, nFire);
lgC     = NaN(nSrc, nDec, nFire);
plgC    = NaN(nSrc, nDec, nFire);
R2      = NaN(nSrc, nDec, nFire);
RMSE    = NaN(nSrc, nDec, nFire);
SSE     = NaN(nSrc, nDec, nFire);
logLik  = NaN(nSrc, nDec, nFire);
AIC     = NaN(nSrc, nDec, nFire);
BIC     = NaN(nSrc, nDec, nFire);

edges = 10 .^ (-1:0.05:5);
for s = 1:nSrc
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{s}));
    sizeAll = [data.size]';

    Dec = {data.decade}'; 
    WUIFire = [data.BrdWUIFire]'; 
    NorCal = [data.NorCal]'; 

    for i = 1:nDec
        for j = 1:nFire
            size = sizeAll(NorCal == 1 & strcmp(Dec, decades{i}) & WUIFire == 2-j);
            N(s,i,j) = numel(size);
    
            if numel(size) < 3, continue; end

            [counts, binEdges] = histcounts(size, edges);
            binWidth = diff(binEdges);
            probDensity = counts ./ (sum(counts) * binWidth); 
            binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
            
            x = log10(binCenters(probDensity > 0)); 
            y = log10(probDensity(probDensity > 0));
            mdl = fitlm(x, y); 

            beta(s,i,j)     = abs(mdl.Coefficients.Estimate(2));
            betaErr(s,i,j)  = mdl.Coefficients.SE(2);
            pBeta(s,i,j)    = mdl.Coefficients.pValue(2);
            lgC(s,i,j)      = mdl.Coefficients.Estimate(1);
            plgC(s,i,j)     = mdl.Coefficients.pValue(1);
            R2(s,i,j)       = mdl.Rsquared.Ordinary;
            RMSE(s,i,j)     = mdl.RMSE;
            SSE(s,i,j)      = mdl.SSE;
            logLik(s,i,j)   = mdl.LogLikelihood;
            AIC(s,i,j)      = mdl.ModelCriterion.AIC;
            BIC(s,i,j)      = mdl.ModelCriterion.BIC;
        end
    end  
end

% Plotting trend line
color = {[0.769, 0.635, 0.820],[0.976, 0.686, 0.545], [0.314, 0.459, 0.314],[0.396, 0.690, 0.910]}; 
Marker = {'d', 'o', '^', 's'};

figure; hold on; set(gcf, 'Position', [100, 100, 550, 400]);
for s = 1:nSrc
    for j = 1:nFire
        y = squeeze(beta(s,:, j));
        e = squeeze(betaErr(s,:, j));
        
        validIdx = ~isnan(y);  
        x = 1:nDec; 
        x = x(validIdx);
        y = y(validIdx);
        e = e(validIdx);

        errorbar(x, y, e, 'Color', color{s}, 'LineStyle', 'none', ...
            'LineWidth', 1.5, 'CapSize', 4, 'HandleVisibility', 'off');

        if j ==1
            plot(x, y, '-o', 'Color', color{s}, 'MarkerFaceColor', color{s}, 'MarkerEdgeColor', 'k', ...
                'Marker', Marker{s}, 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', fires{j});
        elseif j == 2
            plot(x, y, '--o', 'Color', color{s}, 'MarkerFaceColor', color{s}, 'MarkerEdgeColor','k', ...
                'Marker', Marker{s}, 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', fires{j});
        end
    end
end
text(gca, 0.02, 0.05, 'NorCal','Units','normalized', 'HorizontalAlignment','left', ... 
    'VerticalAlignment','bottom', 'FontSize',12, 'FontWeight','bold', ... 
    'BackgroundColor','w', 'Margin',2, 'Clipping','on');
set(gca, 'FontSize', 12);
xlim([0.8, 4.2]); 
ylim([0.95, 2]); 
xticks(1:4); 
xticklabels(decades);
yticks(1:0.2:2);
ylabel('\beta value', 'Interpreter', 'tex', 'FontSize', 14);
set(gca, 'FontSize', 16);
box on;
legend off;

exportgraphics(gcf, 'Figs/Fig3B.pdf', 'ContentType', 'vector');
close(gcf);

% export data: NaN(nSrc, nDec, nFire);
dataSrc_vec = repmat(dataSrc(:), nDec*nFire, 1);                 % 4×(4*2) -> 32×1
decade_vec  = repmat(repelem(decades(:), nSrc, 1), nFire, 1);    % (each decade repeated nSrc) repeated for nFire
fire_vec    = repelem(fires(:), nSrc*nDec, 1);                   % each fire repeated nSrc*nDec
T = table( ...
    dataSrc_vec, decade_vec, fire_vec, ...
    N(:), beta(:), betaErr(:), pBeta(:), lgC(:), plgC(:), R2(:), RMSE(:), SSE(:), logLik(:), AIC(:), BIC(:), ...
    'VariableNames', {'dataSrc','decade','fire','N','Beta','Beta_SE','Beta_p', ...
                      'Log10_C','Log10_C_p','R2','RMSE','SSE','LogLik','AIC','BIC'});
writetable(T, 'Figs/TableS5_Fig3B.csv');
disp(T);
%% SoCal
clear;clc;

dataSrc = {'CalFire', 'MTBS', 'Atlas', 'FIRED'}; nSrc = numel(dataSrc);
decades = {'1990s', '2000s', '2010s', '2020s'}; nDec = numel(decades);
fires = {'WUI', 'Wildland'};                    nFire = numel(fires);

N       = NaN(nSrc, nDec, nFire);
beta    = NaN(nSrc, nDec, nFire);
betaErr = NaN(nSrc, nDec, nFire);
pBeta   = NaN(nSrc, nDec, nFire);
lgC     = NaN(nSrc, nDec, nFire);
plgC    = NaN(nSrc, nDec, nFire);
R2      = NaN(nSrc, nDec, nFire);
RMSE    = NaN(nSrc, nDec, nFire);
SSE     = NaN(nSrc, nDec, nFire);
logLik  = NaN(nSrc, nDec, nFire);
AIC     = NaN(nSrc, nDec, nFire);
BIC     = NaN(nSrc, nDec, nFire);

edges = 10 .^ (-1:0.05:5);
for s = 1:nSrc
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{s}));
    sizeAll = [data.size]';

    Dec = {data.decade}'; 
    WUIFire = [data.BrdWUIFire]'; 
    SoCal = [data.SoCal]'; 

    for i = 1:nDec
        for j = 1:nFire
            size = sizeAll(SoCal == 1 & strcmp(Dec, decades{i}) & WUIFire == 2-j);
            N(s,i,j) = numel(size);
    
            if numel(size) < 3, continue; end

            [counts, binEdges] = histcounts(size, edges);
            binWidth = diff(binEdges);
            probDensity = counts ./ (sum(counts) * binWidth); 
            binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
            
            x = log10(binCenters(probDensity > 0)); 
            y = log10(probDensity(probDensity > 0));
            mdl = fitlm(x, y); 

            beta(s,i,j)     = abs(mdl.Coefficients.Estimate(2));
            betaErr(s,i,j)  = mdl.Coefficients.SE(2);
            pBeta(s,i,j)    = mdl.Coefficients.pValue(2);
            lgC(s,i,j)      = mdl.Coefficients.Estimate(1);
            plgC(s,i,j)     = mdl.Coefficients.pValue(1);
            R2(s,i,j)       = mdl.Rsquared.Ordinary;
            RMSE(s,i,j)     = mdl.RMSE;
            SSE(s,i,j)      = mdl.SSE;
            logLik(s,i,j)   = mdl.LogLikelihood;
            AIC(s,i,j)      = mdl.ModelCriterion.AIC;
            BIC(s,i,j)      = mdl.ModelCriterion.BIC;
        end
    end  
end

% Plotting trend line
color = {[0.769, 0.635, 0.820],[0.976, 0.686, 0.545], [0.314, 0.459, 0.314],[0.396, 0.690, 0.910]}; 
Marker = {'d', 'o', '^', 's'};

figure; hold on; set(gcf, 'Position', [100, 100, 550, 400]);
for s = 1:nSrc
    for j = 1:nFire
        y = squeeze(beta(s,:, j));
        e = squeeze(betaErr(s,:, j));
        
        validIdx = ~isnan(y);  
        x = 1:nDec; 
        x = x(validIdx);
        y = y(validIdx);
        e = e(validIdx);

        errorbar(x, y, e, 'Color', color{s}, 'LineStyle', 'none', ...
            'LineWidth', 1.5, 'CapSize', 4, 'HandleVisibility', 'off');

        if j ==1
            plot(x, y, '-o', 'Color', color{s}, 'MarkerFaceColor', color{s}, 'MarkerEdgeColor', 'k', ...
                'Marker', Marker{s}, 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', fires{j});
        elseif j == 2
            plot(x, y, '--o', 'Color', color{s}, 'MarkerFaceColor', color{s}, 'MarkerEdgeColor','k', ...
                'Marker', Marker{s}, 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', fires{j});
        end
    end
end
text(gca, 0.02, 0.05, 'SoCal','Units','normalized', 'HorizontalAlignment','left', ... 
    'VerticalAlignment','bottom', 'FontSize',12, 'FontWeight','bold', ... 
    'BackgroundColor','w', 'Margin',2, 'Clipping','on');
set(gca, 'FontSize', 12);
xlim([0.8, 4.2]); 
ylim([0.95, 2]); 
xticks(1:4); 
xticklabels(decades);
yticks(1:0.2:2);
ylabel('\beta value', 'Interpreter', 'tex', 'FontSize', 14);
set(gca, 'FontSize', 16);
box on;
legend off;

exportgraphics(gcf, 'Figs/Fig3C.pdf', 'ContentType', 'vector');
close(gcf);

% export data: NaN(nSrc, nDec, nFire);
dataSrc_vec = repmat(dataSrc(:), nDec*nFire, 1);                 % 4×(4*2) -> 32×1
decade_vec  = repmat(repelem(decades(:), nSrc, 1), nFire, 1);    % (each decade repeated nSrc) repeated for nFire
fire_vec    = repelem(fires(:), nSrc*nDec, 1);                   % each fire repeated nSrc*nDec
T = table( ...
    dataSrc_vec, decade_vec, fire_vec, ...
    N(:), beta(:), betaErr(:), pBeta(:), lgC(:), plgC(:), R2(:), RMSE(:), SSE(:), logLik(:), AIC(:), BIC(:), ...
    'VariableNames', {'dataSrc','decade','fire','N','Beta','Beta_SE','Beta_p', ...
                      'Log10_C','Log10_C_p','R2','RMSE','SSE','LogLik','AIC','BIC'});
writetable(T, 'Figs/TableS6_Fig3C.csv');
disp(T);
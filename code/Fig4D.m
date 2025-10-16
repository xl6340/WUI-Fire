%% woodland, grassland, cropland, urban
clear;clc;

groupNames = {'woodland', 'grassland', 'cropland', 'urban'}; nGroup = numel(groupNames);
fires  = {'WUI', 'Wildland'};                                nFire = numel(fires);

N       = NaN(nGroup, nFire); 
beta    = NaN(nGroup, nFire); 
betaErr = NaN(nGroup, nFire); 
pBeta   = NaN(nGroup, nFire); 
lgC     = NaN(nGroup, nFire); 
plgC    = NaN(nGroup, nFire); 
R2      = NaN(nGroup, nFire); 
RMSE    = NaN(nGroup, nFire); 
SSE     = NaN(nGroup, nFire); 
logLik  = NaN(nGroup, nFire); 
AIC     = NaN(nGroup, nFire); 
BIC     = NaN(nGroup, nFire); 

Fcount  = NaN(nGroup, nFire); 
Farea   = NaN(nGroup, nFire); 

edges = 10 .^ (-1:0.05:5); 

data = shaperead(sprintf('dataPrc/firePrmt/CalFire.shp'));
landcover = {data.landcover}';
sizeAll = [data.size]';
WUIFire = [data.BrdWUIFire]'; 

for i = 1:nGroup
    for j = 1:nFire
        sizeGroup = sizeAll(WUIFire == 2-j);
        size = sizeAll(strcmp(landcover, groupNames{i}) & WUIFire == 2-j);
        N(i,j) = numel(size);

        Fcount(i,j) = numel(size) ./ numel(sizeGroup); 
        Farea(i,j)  = sum(size, "all") ./ sum(sizeGroup, "all");

        if numel(size) < 10
            continue;
        end    
        
        [counts, binEdges] = histcounts(size, edges);
        binWidth = diff(binEdges);
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));            
        
        x = log10(binCenters(probDensity > 0)); y = log10(probDensity(probDensity > 0));
        mdl = fitlm(x, y); 

        beta(i,j)    = abs(mdl.Coefficients.Estimate(2));
        betaErr(i,j) = mdl.Coefficients.SE(2);
        pBeta(i,j)   = mdl.Coefficients.pValue(2);
        lgC(i,j)     = mdl.Coefficients.Estimate(1);
        plgC(i,j)    = mdl.Coefficients.pValue(1);
        R2(i,j)      = mdl.Rsquared.Ordinary;
        RMSE(i,j)    = mdl.RMSE;
        SSE(i,j)     = mdl.SSE;
        logLik(i,j)  = mdl.LogLikelihood;
        AIC(i,j)     = mdl.ModelCriterion.AIC;
        BIC(i,j)     = mdl.ModelCriterion.BIC;
    end
end
%% draw bar of beta
colors = [81,130,134; 163,204,208; 224, 238, 172; 203,204,199] / 255;

M      = beta.';     
Merr   = betaErr.';

figure('Position',[100 100 250 300]); hold on
b = bar(M, 'grouped', 'EdgeColor','none', 'BarWidth',1);
for k = 1:numel(b)
    b(k).FaceColor = colors(k,:);
end

for k = 1:numel(b)
    xk = b(k).XEndPoints;
    errorbar(xk, M(:,k), Merr(:,k), 'k', 'linestyle','none', ...
        'LineWidth',1, 'CapSize',4, 'HandleVisibility','off');
end

set(gca,'FontSize',10,'Box','off');
xticks(1:numel(b)); xticklabels({'WUI','Wildland'});
ylim([1 2.5]); yticks(1:0.5:2.5);
ylabel('\beta value');

hL = legend({'Woodland', 'Grassland', 'Cropland', 'Urban'}, 'Location','northwest', 'Box','off');
hL.ItemTokenSize = [10, 25]; 

exportgraphics(gcf, 'Figs/Fig4D1.pdf', 'ContentType', 'vector');
close(gcf);
%% Draw Concentric Cycle (donut) Diagrams
figure('Position',[200 100 160 300]); 
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

r_inner_in  = 0.25;   % inner radius of inner ring (WUI)
r_inner_out = 0.55;   % outer radius of inner ring
r_outer_in  = 0.62;   % inner radius of outer ring (Wildland)
r_outer_out = 1.00;   % outer radius of outer ring

types = {'Fire count','Burned area'};
nan2zeroRow = @(v) (v(:).').*~isnan(v(:).') + 0.*isnan(v(:).');

for t = 1:2
    nexttile; hold on; axis equal off

    switch types{t}
        case 'Fire count'
            vInner = nan2zeroRow(Fcount(:, 1));   % WUI
            vOuter = nan2zeroRow(Fcount(:, 2));   % Wildland
        case 'Burned area'
            vInner = nan2zeroRow(Farea(:, 1));    % WUI
            vOuter = nan2zeroRow(Farea(:, 2));    % Wildland
    end

    draw_annulus_ring(vOuter, colors, r_outer_in, r_outer_out);
    draw_annulus_ring(vInner, colors, r_inner_in, r_inner_out);

    text(0, 1.2, types{t}, 'HorizontalAlignment','center','VerticalAlignment','middle', 'FontWeight','bold');

    y_WUI      = (r_inner_in + r_inner_out)/2;
    y_Wildland = (r_outer_in + r_outer_out)/2;
    text(0.1, y_WUI, 'WUI','HorizontalAlignment','left','VerticalAlignment','top');
    text(0.1, y_Wildland, 'Wildland','HorizontalAlignment','left','VerticalAlignment','top');

    if t == 2
        p = gobjects(1, numel(groupNames));
        for ii = 1:numel(groupNames)
            p(ii) = patch(nan, nan, colors(ii,:), 'EdgeColor','none', 'DisplayName', groupNames{ii});
        end
    end
end

exportgraphics(gcf, 'Figs/Fig4D2.pdf', 'ContentType', 'vector');
close(gcf);

% export data
groups_vec = repmat(groupNames(:), nFire, 1);
fires_vec = repelem(fires(:), nGroup, 1);

T = table( ...
    groups_vec, fires_vec, ...
    N(:), beta(:), betaErr(:), pBeta(:), lgC(:), plgC(:), R2(:), RMSE(:), SSE(:), logLik(:), AIC(:), BIC(:), Fcount(:), Farea(:),...
    'VariableNames', {'groups','fires','N','Beta','Beta_SE','Beta_p', ...
                      'Log10_C','Log10_C_p','R2','RMSE','SSE','LogLik','AIC','BIC', 'Fcount', 'Farea'});
writetable(T, 'Figs/TableS8_Fig4D.csv');
disp(T);
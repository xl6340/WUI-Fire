%% fire spread direction
clear;clc;

fires  = {'WUI', 'Wildland'};                                        nFire = numel(fires);
directions ={'northwest';'southeast';'north';'west';'south';'east'}; nDrctn = numel(directions);

N       = NaN(nDrctn, nFire); 
beta    = NaN(nDrctn, nFire); 
betaErr = NaN(nDrctn, nFire); 
pBeta   = NaN(nDrctn, nFire); 
lgC     = NaN(nDrctn, nFire); 
plgC    = NaN(nDrctn, nFire); 
R2      = NaN(nDrctn, nFire); 
RMSE    = NaN(nDrctn, nFire); 
SSE     = NaN(nDrctn, nFire); 
logLik  = NaN(nDrctn, nFire); 
AIC     = NaN(nDrctn, nFire); 
BIC     = NaN(nDrctn, nFire); 

Fcount   = NaN(nDrctn, nFire); 
Farea    = NaN(nDrctn, nFire); 

data = shaperead('dataPrc/firePrmt/Atlas.shp');
landcover = {data.landcover}';
sizeAll = [data.size]';
WUIFire = [data.BrdWUIFire]'; 
Direction = {data.direction}';

unique(Direction)
edges = 10 .^ (-1:0.05:5); 
for i = 1:nDrctn
    for j = 1:nFire
        sizeGroup = sizeAll(WUIFire == 2-j & ~strcmp(Direction, 'none'));
        size = sizeAll(strcmp(Direction, directions{i}) & WUIFire == 2-j);
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
colors = [ ...
    205 172 180;  % #cdacb4
    233 210 200;  % #e9d2c8
    244 239 203;  % #f4efcb
    214 225 211;  % #d6e1d3
    179 203 226;  % #b3cbe2
    143 171 218   % #8fabda
] / 255;

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
hL = legend(directions, 'Location','northwest', 'Box','off');
hL.ItemTokenSize = [10, 25]; 

exportgraphics(gcf, 'Figs/Fig4B1.pdf', 'ContentType', 'vector');
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
        p = gobjects(1, numel(directions));
        for ii = 1:numel(directions)
            p(ii) = patch(nan, nan, colors(ii,:), 'EdgeColor','none', 'DisplayName', directions{ii});
        end
    end
end

exportgraphics(gcf, 'Figs/Fig4B2.pdf', 'ContentType', 'vector');
close(gcf);


% export data
fires_vec = repelem(fires(:), nDrctn, 1);
direction_vec = repmat(directions(:), nFire, 1);

T = table( ...
    fires_vec, direction_vec, ...
    N(:), beta(:), betaErr(:), pBeta(:), lgC(:), plgC(:), R2(:), RMSE(:), SSE(:), logLik(:), AIC(:), BIC(:), Fcount(:), Farea(:),...
    'VariableNames', {'fires','direction','N','Beta','Beta_SE','Beta_p', ...
                      'Log10_C','Log10_C_p','R2','RMSE','SSE','LogLik','AIC','BIC', 'Fcount', 'Farea'});
writetable(T, 'Figs/TableS6_Fig4B.csv');
disp(T);


%% ignition source, CalFire
clear;clc;

Nature = [1 17]; Human = 2:19;
groupNames = {'Nature', 'Human'};
groups  = {Nature, Human};         nGroup = numel(groupNames);
fires   = {'WUI', 'Wildland'};     nFire = numel(fires);
decades = {'Pre2010', 'Post2010'}; nDec = numel(decades);

N       = NaN(nGroup, nFire, nDec); 
beta    = NaN(nGroup, nFire, nDec); 
betaErr = NaN(nGroup, nFire, nDec); 
pBeta   = NaN(nGroup, nFire, nDec); 
lgC     = NaN(nGroup, nFire, nDec); 
plgC    = NaN(nGroup, nFire, nDec); 
R2      = NaN(nGroup, nFire, nDec); 
RMSE    = NaN(nGroup, nFire, nDec); 
SSE     = NaN(nGroup, nFire, nDec); 
logLik  = NaN(nGroup, nFire, nDec); 
AIC     = NaN(nGroup, nFire, nDec); 
BIC     = NaN(nGroup, nFire, nDec); 

Fcount  = NaN(nGroup, nFire, nDec); 
Farea   = NaN(nGroup, nFire, nDec); 

data = shaperead(sprintf('dataPrc/firePrmt/CalFire.shp'));
cause = [data.CAUSE]';
sizeAll = [data.size]';
WUIFire = [data.BrdWUIFire]';
decade = {data.decade}'; 
edges = 10 .^ (-1:0.05:5); 

for i = 1:nGroup
    for j = 1:nFire
        for d = 1:nDec
            if d == 1, mask = ismember(decade, {'1990s', '2000s'});
            else, mask = ismember(decade, {'2010s', '2020s'}); end

            sizeGroup = sizeAll(WUIFire == 2-j & mask);
            size = sizeAll(ismember(cause, groups{i}) & WUIFire == 2-j & mask);
            N(i,j,d) = numel(size);
    
            Fcount(i,j,d) = numel(size) ./ numel(sizeGroup); 
            Farea(i,j,d)  = sum(size, "all") ./ sum(sizeGroup, "all");
    
            if numel(size) < 10
                continue;
            end    
            
            [counts, binEdges] = histcounts(size, edges);
            binWidth = diff(binEdges);
            probDensity = counts ./ (sum(counts) * binWidth); 
            binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));            
            
            x = log10(binCenters(probDensity > 0)); y = log10(probDensity(probDensity > 0));
            mdl = fitlm(x, y); 
    
            beta(i,j,d)    = abs(mdl.Coefficients.Estimate(2));
            betaErr(i,j,d) = mdl.Coefficients.SE(2);
            pBeta(i,j,d)   = mdl.Coefficients.pValue(2);
            lgC(i,j,d)     = mdl.Coefficients.Estimate(1);
            plgC(i,j,d)    = mdl.Coefficients.pValue(1);
            R2(i,j,d)      = mdl.Rsquared.Ordinary;
            RMSE(i,j,d)    = mdl.RMSE;
            SSE(i,j,d)     = mdl.SSE;
            logLik(i,j,d)  = mdl.LogLikelihood;
            AIC(i,j,d)     = mdl.ModelCriterion.AIC;
            BIC(i,j,d)     = mdl.ModelCriterion.BIC;
        end
    end
end
%% draw bar of beta (two panels side by side)
colors = [150 210 176; 243 215 138] ./ 255; 

M1    = beta(:,:,1).';     
Merr1 = betaErr(:,:,1).';
M2    = beta(:,:,2).';     
Merr2 = betaErr(:,:,2).';

M_all    = [M1; M2]; % Combine both groups horizontally
Merr_all = [Merr1; Merr2];

figure('Position',[100 100 250 300]); hold on

b = bar(M_all, 'grouped', 'EdgeColor','none', 'BarWidth', 1);
for k = 1:numel(b)
    b(k).FaceColor = colors(k,:);
end

for k = 1:numel(b)
    xk = b(k).XEndPoints;
    errorbar(xk, M_all(:,k), Merr_all(:,k), 'k', 'linestyle','none', ...
        'LineWidth', 1, 'CapSize', 4, 'HandleVisibility', 'off');
end

xline(2.5, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);

set(gca,'FontSize',10,'Box','off');
xticks(1:4); 
xticklabels({'WUI','Wildland','WUI','Wildland'});
xtickangle(30);
ylim([0.8 2]); yticks(0.8:0.2:2);
ylabel('\beta value');

hL = legend({'Nature','Human'}, 'Location','best', 'Box','off');
hL.ItemTokenSize = [10, 25];
pos = hL.Position; pos(2) = pos(2) - 0.08; hL.Position = pos;

yl = ylim;
text(1.5, yl(2) - 0.05, 'pre2010', ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight','bold');
text(3.5, yl(2) - 0.05, 'post2010', ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight','bold');

exportgraphics(gcf, 'Figs/Fig4C1.pdf', 'ContentType', 'vector');
close(gcf);
%% Draw Concentric Cycle (donut) Diagrams
colNature = [150 210 176] / 255;   % soft green
colHuman  = [243 215 138] / 255;   % soft yellow
sliceColors = [colNature; colHuman];

figure('Position',[200 100 300 160]); 
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

Radii for the rings
r_inner_in  = 0.25;   % inner radius of inner ring (WUI)
r_inner_out = 0.55;   % outer radius of inner ring (WUI)
r_outer_in  = 0.62;   % inner radius of outer ring (Wildland)
r_outer_out = 1.00;   % outer radius of outer ring (Wildland)

Helpers
types = {'Count','Area'};
nan2zero = @(v) (v.' .* ~isnan(v.')) + 0 .* isnan(v.');

for t = 1:2
    nexttile; hold on; axis equal off
    switch types{t}
        case 'Count'
            vInner = nan2zero(Fcount(:, 1));   % WUI  (inner ring)
            vOuter = nan2zero(Fcount(:, 2));   % Wildland (outer ring)
        case 'Area'
            vInner = nan2zero(Farea(:,  1));   % WUI  (inner ring)
            vOuter = nan2zero(Farea(:,  2));   % Wildland (outer ring)
    end

    draw_annulus_ring(vOuter, sliceColors, r_outer_in, r_outer_out);
    draw_annulus_ring(vInner, sliceColors, r_inner_in, r_inner_out);

    text(0, 0, types{t}, 'HorizontalAlignment','center', ...
        'VerticalAlignment','middle','FontSize',9, 'Color','k');

    y_WUI      = - (r_inner_in + r_inner_out)/2;   % inner ring midpoint
    y_Wildland = - (r_outer_in + r_outer_out)/2;   % outer ring midpoint

    text(0, y_WUI, 'WUI', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',9, 'Color','k');
    text(0, y_Wildland, 'Wildland', ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',9, 'Color','k');

    if t == 2
        patch(nan, nan, colNature, 'EdgeColor','none', 'DisplayName','Nature');
        patch(nan, nan, colHuman,  'EdgeColor','none', 'DisplayName','Human');
    end
end

% exportgraphics(gcf, 'Figs/Fig4C2.pdf', 'ContentType', 'vector');
% close(gcf);


% export data
groups_vec = repelem(groupNames(:), nFire, 1);
fires_vec = repmat(fires(:), nGroup, 1);

for t=1:2
     N1       = N(:,:,t); 
     beta1    = beta(:,:,t); 
     betaErr1 = betaErr(:,:,t); 
     pBeta1   = pBeta(:,:,t); 
     lgC1     = lgC(:,:,t); 
     plgC1    = plgC(:,:,t); 
     R21      = R2(:,:,t); 
     RMSE1    = RMSE(:,:,t); 
     SSE1     = SSE(:,:,t); 
     logLik1  = logLik(:,:,t); 
     AIC1     = AIC(:,:,t); 
     BIC1     = BIC(:,:,t); 
     Fcount1  = Fcount(:,:,t); 
     Farea1   = Farea(:,:,t); 
    
    T = table(groups_vec, fires_vec, ...
        N1(:), beta1(:), betaErr1(:), pBeta1(:), lgC1(:), plgC1(:), R21(:), RMSE1(:), SSE1(:), logLik1(:), AIC1(:), BIC1(:), Fcount1(:), Farea1(:),...
        'VariableNames', {'groups','fires','N','Beta','Beta_SE','Beta_p', ...
                          'Log10_C','Log10_C_p','R2','RMSE','SSE','LogLik','AIC','BIC', 'Fcount', 'Farea'});
    writetable(T, sprintf('Figs/TableS7_Fig4C_%s.csv',decades{t}));
end


Nature = [1 17]; Human = 2:19;
groupNames = {'Nature', 'Human'};
groups  = {Nature, Human};         nGroup = numel(groupNames);
fires   = {'WUI', 'Wildland'};     nFire = numel(fires);
decades = {'Pre2010', 'Post2010'}; nDec = numel(decades);

N       = NaN(nGroup, nFire, nDec); 
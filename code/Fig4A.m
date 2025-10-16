%% climate, CalFire
clear; clc;

prisms = {'ppt30sum','vpd3mean'}; nPrism = numel(prisms);
fires   = {'WUI', 'Wildland'};    nFire = numel(fires);

data       = shaperead('/Users/xinlei/Desktop/firePerimeter/CalFire.shp');

WUIFire = [data.BrdWUIFire]'; 
pptAll = [data.ppt30sum]'; 
vpdAll    = [data.vpd3mean]'; 

% axis labels
xlabels.ppt    = 'PPT (log_{10}, mm)';
xlabels.vpdmax = 'VPD_{max} (hPa)';
xlab_list = {xlabels.ppt, xlabels.vpdmax}; 

figure; clf; set(gcf,'Position',[100, 100, 400, 300]);
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

colr = {[0.8 0.4 0.2], [0.2 0.6 0.8]};  % {WUI, Wildland}

for i = 1:nPrism
    if i == 1
        prismWUI  = log10(pptAll(WUIFire == 1 & isfinite(pptAll) & pptAll > 0));
        prismWild = log10(pptAll(WUIFire == 0 & isfinite(pptAll) & pptAll > 0));
    elseif i==2
        prismWUI  = vpdAll(WUIFire == 1 & isfinite(vpdAll) & vpdAll > 0);
        prismWild = vpdAll(WUIFire == 0 & isfinite(vpdAll) & vpdAll > 0);
    end

    [~, p]   = ttest2(prismWUI, prismWild);
    meanWUI  = mean(prismWUI,'omitnan');  stdWUI  = std(prismWUI,'omitnan');
    meanWild = mean(prismWild,'omitnan'); stdWild = std(prismWild,'omitnan');

    ax = nexttile(i); hold(ax,'on');
    
    edges = linspace(min([prismWUI(:); prismWild(:)]), max([prismWUI(:); prismWild(:)]), 20);
    histogram(prismWUI, 'BinEdges', edges, 'Normalization', 'probability', ...
        'FaceColor', colr{1}, 'FaceAlpha', 0.5, 'EdgeColor', colr{1});
    hold on
    histogram(prismWild, 'BinEdges', edges, 'Normalization', 'probability', ...
        'FaceColor', colr{2}, 'FaceAlpha', 0.5, 'EdgeColor', colr{2});

    % histogram(prismWUI,'Normalization','probability',...
    %     'FaceColor',colr{1},'FaceAlpha',0.5,'EdgeColor',colr{1});
    % histogram(prismWild,'Normalization','probability',...
    %     'FaceColor',colr{2},'FaceAlpha',0.5,'EdgeColor',colr{2});

    yl = ylim; xl = xlim;
    plot([meanWUI meanWUI], yl, '--', 'Color', colr{1},'HandleVisibility','off', 'LineWidth',1.5);
    plot([meanWild meanWild], yl, '--', 'Color', colr{2},'HandleVisibility','off', 'LineWidth',1.5);

    if p < 0.01, ptxt = '\it{p} < 0.01';
    else,        ptxt = sprintf('\\it{p} = %.3f', p);
    end
    text(xl(2)*0.98, yl(2)*0.95, {ptxt}, 'HorizontalAlignment','right','VerticalAlignment','top','Interpreter','tex');
    xlabel(xlab_list{i});
    ylabel('Probability');
    lgd = legend({'WUI','Wildland'}, 'Location','northwest', 'Box','off');
    lgd.ItemTokenSize(1) = 9;
    set(ax,'FontSize',9,'Box','off', 'LineWidth', 0.25);
    hold off;
end

exportgraphics(gcf, 'Figs/Fig4A.pdf', 'ContentType', 'vector');
close(gcf);

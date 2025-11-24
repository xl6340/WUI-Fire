%% Figure 4. NDVI
clear; clc;

color = {[216, 118, 89]/255; [41, 157, 143]/255}; % Urban, Wildland
FireType = {'Urban-edge', 'Wildland'};  nFire = numel(FireType);
decades  = {'2000s', '2010s', '2020s'}; nDec = numel(decades);

figure('Position',[100 100 400 600]);
t = tiledlayout(nDec, 1, 'TileSpacing', 'compact');
axes_handles = gobjects(nDec, 1);
for d = 1:nDec
    ax = nexttile;
    axes_handles(d) = ax;
    hold(ax,'on');
    set(ax,'TickDir', 'out', 'FontSize', 10);

    ndviU = []; ndviW = [];
    edge = 0:0.05:1;
    for j = nFire:-1:1  
        data = readtable(sprintf('dataFig/ndvi/%s-%s.csv', FireType{j}, decades{d}));
        h = histogram(ax, [data.ndvi], 'BinEdges', edge,'Normalization','probability', 'FaceColor', color{j}, ...
            'EdgeColor', [1 1 1], 'FaceAlpha', 0.4,'DisplayName', FireType{j});

        xline(ax, nanmean([data.ndvi]), 'LineWidth', 2.5, 'LineStyle', ':', ...
                'Color', color{j}, 'HandleVisibility', 'off'); 

        if j == 1,ndviU = [data.ndvi]; else,ndviW = [data.ndvi]; end
    end
    ax.TickDir = 'out';
    xlim(ax, [0 1]); 
    xticks(ax, 0:0.2:1);
    ylim(ax, [0 0.2]);
    yticks(ax, 0:0.1:0.2);

    % if d == 1 
    %     lg = legend(ax, 'show','Location', 'southeast', 'FontSize',14);
    %     % lg.ItemTokenSize = [12 10]; 
    % end

    yl = ylim(ax); 
    xl = xlim(ax);

    text(ax, xl(2) - 0.02*diff(xl), yl(2)*0.90, decades{d}, ...
        'Color', [0 0 0], 'FontSize', 12, 'HorizontalAlignment', 'right');

    text(ax, xl(1)+0.02*diff(xl), yl(2)*0.90, ...
        sprintf('Mean: %.2f \\pm %.2f', nanmean(ndviU), nanstd(ndviU)), ...
        'Color', color{1}, 'FontSize',10);
    text(ax, xl(1)+0.02*diff(xl), yl(2)*0.75, ...
        sprintf('Mean: %.2f \\pm %.2f', nanmean(ndviW), nanstd(ndviW)), ...
        'Color', color{2}, 'FontSize',10);
    
    [~, p_ttest] = ttest2(ndviU, ndviW, 'Vartype','unequal');
    text(ax, xl(1)+0.02*diff(xl), yl(2)*0.60, sprintf('\\it{p}\\rm{} = %.4f', p_ttest), ...
        'Color', [0 0 0], 'Interpreter','tex', 'FontSize',10);
end
xlabel(t, 'NDVI', 'Interpreter', 'latex', 'FontSize',14);
ylabel(t, 'Probability', 'FontSize',14);
linkaxes(axes_handles, 'xy');

exportgraphics(gcf,'dataFig/Figs/Fig4.pdf','ContentType','vector');
close(gcf);
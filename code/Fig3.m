%% Figure 3. climate
clear; clc;
varNames      = {'tmean', 'vpdmax','tmin', 'tmax'};
xLabels       = {'Tmean (°C)', 'VPDmax (hPa)','Tmin (°C)', 'Tmax (°C)'};
xLimits       = {[0 40], [0 80], [-5 35], [5 45]}; 
xTicks        = {-5:5:40, 0:10:80, -5:5:35, 5:5:45}; 
yLimits       = {[0 0.15], [0 0.1], [0 0.15], [0 0.15]}; 
yTicks        = {0:0.05:0.15, 0:0.02:0.1, 0:0.05:0.15,0:0.05:0.15}; 
binWidths     = {1.2, 2, 1.2, 1.2};
units         = {'°C', 'hPa', '°C', '°C'};
nVars = numel(varNames);
color = {[216, 118, 89]/255; [41, 157, 143]/255}; % Urban, Wildland
FireType = {'Urban-edge', 'Wildland'}; 
nFire = numel(FireType);

data_Urban = readtable('dataFig/climate/Urban-edge.csv');
data_Wild = readtable('dataFig/climate/Wildland.csv');

figure('Position',[100 100 800 400]);
labels = {'(A)', '(B)', '(C)', '(D)'}; 
t = tiledlayout(2, 2);
for i = 1:nVars
    ax = nexttile; hold(ax,'on');  

    dataU = data_Urban.(varNames{i}); dataU = dataU(~isnan(dataU));
    histogram(ax, dataU, 'Normalization','probability', ...
        'FaceColor', color{1}, 'EdgeColor',[1 1 1], 'LineWidth', 0.1, ...
        'BinWidth', binWidths{i},'DisplayName', FireType{1}, 'FaceAlpha', 0.4); 

    dataW = data_Wild.(varNames{i}); dataW = dataW(~isnan(dataW));
    histogram(ax, dataW, 'Normalization','probability', ...
        'FaceColor', color{2}, 'EdgeColor',[1 1 1], 'LineWidth', 0.1,...
        'BinWidth', binWidths{i},'DisplayName', FireType{2},'FaceAlpha', 0.4); 

    xline(ax, nanmean(dataU), 'LineWidth', 1.5, 'LineStyle', '--', ...
        'Color', color{1}, 'HandleVisibility', 'off');
    xline(ax, nanmean(dataW), 'LineWidth', 1.5, 'LineStyle', '--', ...
        'Color', color{2}, 'HandleVisibility', 'off');

    ylim(ax, yLimits{i}); 
    yticks(ax, yTicks{i});
    
    ax.TickDir = 'out';
    
    xlim(ax, xLimits{i}); 
    xticks(ax, xTicks{i});
    ylabel(ax, 'Probability', 'FontSize', 12);
    xlabel(ax, xLabels{i}, 'FontSize', 12);    
    set(ax, 'FontSize', 10);
    
    if i == 1 
        lg = legend(ax, 'show', 'Box', 'off','Location', 'northeast');
        lg.ItemTokenSize = [10 10]; 
    end
    yl = ylim(ax); xl = xlim(ax);
    text(xl(1)+0.02*diff(xl), yl(2)*0.90, ...
        sprintf('Mean: %.1f ± %.1f %s', nanmean(dataU), ...
        nanstd(dataU), units{i}),'Color', color{1}, 'FontSize',10);
    text(xl(1)+0.02*diff(xl), yl(2)*0.75, sprintf('Mean: %.1f ± %.1f %s', ...
        nanmean(dataW), nanstd(dataW), units{i}),'Color', color{2}, 'FontSize',10);
    [~, p_ttest] = ttest2(dataU, dataW, 'Vartype','unequal');
    text(xl(1)+0.02*diff(xl), yl(2)*0.60, sprintf('\\it{p}\\rm{} = %.4f', p_ttest), ...
        'Color', [0 0 0], 'Interpreter','tex', 'FontSize',10);  
    hold(ax, 'off');

    pos = get(ax, 'Position');
    label_x_start = pos(1) - 0.02; 
    label_y_start = pos(2) + pos(4) + 0.01;   
    annotation('textbox', [label_x_start, label_y_start, 0, 0], ...
               'String', labels{i}, 'EdgeColor', 'none', ...
               'FontSize',9, 'HorizontalAlignment', 'left', ...
               'VerticalAlignment', 'bottom', 'FitBoxToText', 'on');
end

exportgraphics(gcf,'dataFig/Figs/Fig3.pdf','ContentType','vector');
close(gcf);
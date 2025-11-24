clear; clc;

dataSrc  = {'CalFire', 'MTBS', 'Atlas', 'FIRED'};   nSrc  = numel(dataSrc);
fireType = {'Urban-edge','Wildland'};               nFire = numel(fireType);
igType   = {'Human', 'Natural'};                    nIg   = numel(igType);
decades  = {'1990s', '2000s', '2010s', '2020s'};    nDec  = numel(decades);
years    = 1990:2024;                               nyr   = numel(years);
lable    = {'A', 'B', 'C', 'D', 'E', 'F'};          nLbl  = numel(lable);

colors = {[216, 118, 89;242, 193, 78]/255; % urban: human, natural
    [41, 157, 143; 138, 201, 38]/255       % wildland: human, natural
    };

figure('Position',[100 100 700 700]);

% Subplot 1: fire number, Urban-edge % [left, bottom, width, height]
ax1 = axes('Position', [0.08, 0.7, 0.38, 0.2]); hold on;
p1 = get(ax1, 'Position'); 
data = readtable(sprintf('dataFig/ba/count-%s.csv', fireType{1}));  
x = data.Row; xfit = linspace(min(x), max(x), 100);
for i = 1:nIg      
    y = data.(igType{i});
    x_design = [ones(size(x)) x];
    [b,~,~,~,stats] = regress(y, x_design);
    yfit = b(1) + b(2) * xfit;

    plot(x, y, '-o', 'Color', [colors{1}(i,:), 0.5], 'MarkerFaceColor', colors{1}(i,:), 'MarkerSize',5);    
    
    if stats(3) <  0.01
        plot(xfit, yfit, '-', 'Color', [0.7 0.7 0.7]);
        text(0.1, 0.5, sprintf('$s = %.2f, p = %.3f$', b(2), stats(3)),...
            'Interpreter', 'latex', 'Units','normalized','Color',colors{1}(i,:),'FontSize',9);
    else
        plot(xfit, yfit, '-.', 'Color', [0.7 0.7 0.7]);        
    end
end
ylabel('Fire number');ylim([0 400]); 
xlim([1990 2025]); xticks(1990:5:2025);xticklabels([]);
set(gca,'FontSize',9); 
title('Urban-edge'); box on;hold off;

% Subplot 2: fire number, Wildland
ax2 = axes('Position', [0.55, 0.7, 0.38, 0.2]); hold on;
p2 = get(ax2, 'Position'); 
data = readtable(sprintf('dataFig/ba/count-%s.csv', fireType{2}));  
x = data.Row; xfit = linspace(min(x), max(x), 100);
for i = 1:nIg      
    y = data.(igType{i});
    x_design = [ones(size(x)) x];
    [b,~,~,~,stats] = regress(y, x_design);
    yfit = b(1) + b(2) * xfit;

    plot(x, y, '-o', 'Color', [colors{2}(i,:), 0.5], 'MarkerFaceColor', colors{2}(i,:), 'MarkerSize',5);
    
    
    if stats(3) <  0.01
        plot(xfit, yfit, '-', 'Color', [0.7 0.7 0.7]);
        text(0.1, 0.6, sprintf('$s = %.2f, p = %.3f$', b(2), stats(3)),...
            'Interpreter', 'latex', 'Units','normalized','Color',colors{2}(i,:),'FontSize',9);
    else
        plot(xfit, yfit, '-.', 'Color', [0.7 0.7 0.7]);
    end
end
ylim([0 400]); 
xlim([1990 2025]);xticks(1990:5:2025);xticklabels([]);
set(gca,'FontSize',9); box on;
title('Wildland'); hold off;

% Subplot 3: burned area, Urban-edge % [left, bottom, width, height]
ax3 = axes('Position', [0.08, 0.42, 0.38, 0.2]); hold on;
p3 = get(ax3, 'Position'); 
data = readtable(sprintf('dataFig/ba/area-%s.csv', fireType{1}));  
x = data.Row; xfit = linspace(min(x), max(x), 100);
for i = 1:nIg      
    y = data.(igType{i});
    x_design = [ones(size(x)) x];
    [b,~,~,~,stats] = regress(y, x_design);
    yfit = b(1) + b(2) * xfit;

    plot(x, y, '-o', 'Color', [colors{1}(i,:), 0.5], 'MarkerFaceColor', colors{1}(i,:), 'MarkerSize',5);    
    
    if stats(3) <  0.01
        plot(xfit, yfit, '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
        text(0.1, 0.4, sprintf('$s = %.2f, p = %.3f$', b(2), stats(3)),...
            'Interpreter', 'latex', 'Units','normalized','Color',colors{1}(i,:),'FontSize',9);
    else
        plot(xfit, yfit, '-.', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    end
end
lgd = legend('Human ignited', 'Natural ignited','box', 'off',...
    'Location', 'north', 'FontSize',9);
lgd.ItemTokenSize = [12 10]; 

xlim([1990 2025]); xticks(1990:5:2025); xlabel('Year');
ylabel('Burned area (×10^3 km^2)'); 
ylim([0 12000]);yticks(0:4000:12000);
yticklabels(arrayfun(@(x) sprintf('%d', x/1000), yticks, 'UniformOutput', false));
set(gca,'FontSize',9); box on;
hold off;

% Subplot 4: burned area, Wildland
ax4 = axes('Position', [0.55, 0.42, 0.38, 0.2]); hold(ax4, 'on');
p4 = get(ax4, 'Position'); 
data = readtable(sprintf('dataFig/ba/area-%s.csv', fireType{2}));  
x = data.Row; xfit = linspace(min(x), max(x), 100);
for i = 1:nIg      
    y = data.(igType{i});
    x_design = [ones(size(x)) x];
    [b,~,~,~,stats] = regress(y, x_design);
    yfit = b(1) + b(2) * xfit;

    plot(ax4, x, y, '-o', 'Color', [colors{2}(i,:), 0.5], 'MarkerFaceColor', colors{2}(i,:), 'MarkerSize',5);    
    
    if stats(3) <  0.01
        plot(ax4, xfit, yfit, '-', 'Color', [0.7 0.7 0.7],'HandleVisibility','off');
        text(0.1, 0.5, sprintf('$s = %.2f, p = %.3f$', b(2), stats(3)),...
            'Interpreter', 'latex', 'Units','normalized','Color',colors{2}(i,:),'FontSize',9);
    else
        plot(ax4, xfit, yfit, '-.', 'Color', [0.7 0.7 0.7],'HandleVisibility','off');
    end
end
lgd = legend(ax4,'Human ignited', 'Natural ignited', 'box', 'off', 'Location', 'north', 'FontSize',9);
lgd.ItemTokenSize = [12 10]; 

xlim([1990 2025]); xticks(1990:5:2025); xlabel('Year');
ylim([0 12000]);yticks(0:4000:12000);
yticklabels(arrayfun(@(x) sprintf('%d', x/1000), yticks, 'UniformOutput', false));
set(gca,'FontSize',9); box on;
hold off;

% Subplot 5: beta of 4 datasets
ax5 = axes('Position', [0.08, 0.1, 0.27, 0.2]); hold on; 
p5 = get(ax5, 'Position'); 
barWidth = 0.25;
for i = 1:nFire    
    data = readtable(sprintf('dataFig/beta/4Srcs-%s.csv', fireType{i}));  
    x = 1:nSrc; y = [data.beta]'; e = [data.betaErr]';
    if i == 1, bar_x = x - barWidth / 2;  
    else, bar_x = x + barWidth / 2; end
    bar(bar_x, y, 'FaceColor', colors{1}(i,:), 'BaseValue', -2, 'BarWidth', barWidth);        
    errorbar(bar_x, y, e, 'Color', 'k', 'LineStyle', ...
        'none', 'CapSize', 4, 'HandleVisibility', 'off');
end
set(gca,'FontSize',9); box on;
xlim([0.3, 4.7]); xticklabels(dataSrc);xtickangle(25); 
ylim([-2, -1]);yticks(-2:0.2:-1); ylabel('\beta value', 'Interpreter', 'tex');
h = legend('Urban-edge', 'Wildland', 'box', 'off', 'Location','northwest', 'FontSize',9);
h.ItemTokenSize = [10, 20];

% Subplot 6: beta of 4 decades, % [left, bottom, width, height]
ax6 = axes('Position', [0.43, 0.1, 0.27, 0.2]); hold on; 
p6 = get(ax6, 'Position'); hold on; 
edgeCut  = {'1-10 km^2', '1-100 km^2', '1-1000 km^2', 'all >1 km^2'}; 
nEdge    = numel(edgeCut);
colorCut = {
    [232, 136, 101;  96, 184, 171] / 255;   % cut 1
    [236, 190, 178;  205, 224, 199] / 255;   % cut 2
    [227, 164, 151;  155, 184, 156] / 255;   % set 3
    [216 118  89;  41 157 143] / 255;   % cut 4
}; % Urban, Wildland 
Marker = {'o', '^'};
plotHandles = gobjects(nEdge, nFire);
for i = 2:nEdge
    for j = 1:nFire
        data = readtable(sprintf('dataFig/beta/4Decs-%s-%s.csv', edgeCut{i}, fireType{j})); 
        x = 1:nDec; y = [data.beta]'; e = [data.betaErr]';
        h_plot = plot(x, y, 'Color', colorCut{i}(j,:), 'MarkerFaceColor', colorCut{i}(j,:), 'MarkerEdgeColor', colorCut{i}(j,:)*0.9, ...
            'Marker', Marker{j},'LineWidth', 1,'MarkerSize', 5);
        errorbar(x, y, e, 'Color', colorCut{i}(j,:), 'LineStyle', 'none', ...
            'LineWidth', 1.5, 'CapSize', 4, 'HandleVisibility', 'off');
        plotHandles(i, j) = h_plot;
    end
end
set(gca,'FontSize',9); box on;
xlim([0.8, 4.2]); xticks(1:4); 
xticklabels(decades);
xtickangle(25); 
ylim([-2, -1]); yticks(-2:0.2:-1);
% Urban legend
urbanHandles = plotHandles(nEdge:-1:2, 1);
urbanLabels  = edgeCut(nEdge:-1:2);
emptyLabels  = repmat({''}, size(urbanLabels)); 
lgd1 = legend(urbanHandles, emptyLabels, ...
    'Location','southeast', 'Box','off', 'FontSize',9);
lgd1.ItemTokenSize = [12 8];
% Wildland legend
ax_top = axes('Position', get(gca,'Position'), 'Visible','off');
wildlandHandles = plotHandles(nEdge:-1:2, 2);
% wildlandLabels  = edgeCut(nEdge:-1:2);
wildlandLabels = {'original', '<1000 km^2 fires', '<100 km^2 fires'};
lgd2 = legend(ax_top, wildlandHandles, wildlandLabels, ...
    'Location','southeast', 'Box','off','FontSize',9);
lgd2.ItemTokenSize = [12 8];

drawnow;
lgd1.Units = 'normalized'; lgd2.Units = 'normalized';
p1 = lgd1.Position; p2 = lgd2.Position;
y_shift_down = 0.01;
p2(2) = p2(2) - y_shift_down; lgd2.Position = p2;
p1 = lgd1.Position; shift = 0.03; 
lgd1.Position = [p2(1) - shift, p2(2),p1(3), p2(4)];


% Subplot 7: beta of 2 ignitions, % [left, bottom, width, height]
ax7 = axes('Position', [0.78, 0.10, 0.15, 0.2]); hold on; 
p7 = get(ax7, 'Position'); hold on; 
for i = 1:nIg   
    data = readtable(sprintf('dataFig/beta/2Fires-%s.csv', igType{i})); 
    x = 1:nIg; y = [data.beta]'; e = [data.betaErr]';    
    if i == 1, bar_x = x - barWidth / 2;  
    else, bar_x = x + barWidth / 2; end
    h = bar(bar_x, y, 'BaseValue', -2, 'BarWidth', barWidth);    
    h.FaceColor = 'flat';
    h.CData = [
        colors{1}(i,:);   % bar 1: urban color for this ignition
        colors{2}(i,:)    % bar 2: wildland color for this ignition
    ];
    for k = 1:numel(y) 
        errorbar(bar_x(k), y(k), e(k),'LineStyle', 'none', ...
            'CapSize', 4, 'Color', 'k');
    end
end
set(gca,'FontSize',9); box on;
xticks(1:nFire); xlim([0.3, 2.7]); 
xticklabels(fireType); xtickangle(25); 
ylim([-2, -1]);yticks(-2:0.2:-1); 

% labels
labels = {'(A)', '(B)', '(C)', '(D)', '(E)', '(F)', '(G)'}; % Assuming A-G for the 7 plots
axes_list = {ax1, ax2, ax3, ax4, ax5, ax6, ax7};

for k = 1:numel(axes_list)
    ax = axes_list{k};
    pos = get(ax, 'Position');    
    label_x_start = pos(1) - 0.02;     
    label_y_start = pos(2) + pos(4) + 0.01;     
    annotation('textbox', [label_x_start, label_y_start, 0, 0], ...
               'String', labels{k}, 'EdgeColor', 'none', ...
               'FontSize',9, 'HorizontalAlignment', 'left', ...
               'VerticalAlignment', 'bottom', 'FitBoxToText', 'on');
end

exportgraphics(gcf,'dataFig/Figs/Fig2.pdf','ContentType','vector');
close(gcf);
%% beta for 4-datasets, 4-decades
clear; clc;

dataSrc  = {'CalFire', 'MTBS', 'Atlas', 'FIRED'};   nSrc  = numel(dataSrc);
decades  = {'1990s', '2000s', '2010s', '2020s'};    nDec  = numel(decades);

count    = NaN(nSrc, nDec);
area     = NaN(nSrc, nDec);
alfa     = NaN(nSrc, nDec);
beta     = NaN(nSrc, nDec);
betaErr  = NaN(nSrc, nDec);
pBeta    = NaN(nSrc, nDec);
R2       = NaN(nSrc, nDec);

edges = 10 .^ (-1:0.05:5); 
for i = 1:nSrc
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{i}));
    Size = [data.size]';
    Dec  = {data.decade}';
    for d = 1:nDec
        sz = Size(strcmp(Dec, decades{d}));
        if isempty(sz), continue; end
        [counts, binEdges] = histcounts(sz, edges);
        binWidth = diff(binEdges);
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
    
        x = log10(binCenters(probDensity > 0)); 
        y = log10(probDensity(probDensity > 0));
        writetable(table(binCenters', probDensity', 'VariableNames', {'binCenter', 'probDensity'}), ...
            sprintf('dataFig/curve/%s-%s.csv', dataSrc{i}, decades{d}));

        mdl = fitlm(x, y);
        x_fit = linspace(min(x), max(x), 100)';
        y_fit = predict(mdl, x_fit);
        writetable(table(10.^(x_fit), 10.^(y_fit), 'VariableNames', {'x_fit', 'y_fit'}), ...
            sprintf('dataFig/curve/%s-%s-fit.csv', dataSrc{i}, decades{d}));

        count(i, d)   = numel(sz);
        area(i, d)    = sum(sz);
        alfa(i, d)    = mdl.Coefficients.Estimate(1);
        beta(i, d)    = mdl.Coefficients.Estimate(2);
        betaErr(i, d) = mdl.Coefficients.SE(2);
        pBeta(i, d)   = mdl.Coefficients.pValue(2);
        R2(i, d)      = mdl.Rsquared.Ordinary;
    end
end

for i = 1:nSrc
        countVec = squeeze(count(i, :))';
        areaVec = squeeze(area(i, :))';
        alfaVec = squeeze(alfa(i, :))';
        betaVec = squeeze(beta(i, :))'; 
        betaErrVec = squeeze(betaErr(i, :))';
        pBetaVec = squeeze(pBeta(i, :))';
        R2Vec = squeeze(R2(i, :))';

        writetable(table(countVec, areaVec, alfaVec, betaVec, betaErrVec, pBetaVec, R2Vec, ...
            'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2'}, ...
            'RowNames', decades), sprintf('dataFig/beta/4Decs-%s.csv', dataSrc{i}), 'WriteRowNames', true);
end
%% beta for 4-datasets,2-fireTypes
clear; clc;

dataSrc  = {'CalFire', 'MTBS', 'Atlas', 'FIRED'};   nSrc  = numel(dataSrc);
fireType = {'Urban-edge','Wildland'};               nFire = numel(fireType);

edges = 10 .^ (-1:0.05:5); 
count    = NaN(nSrc, nFire);
area     = NaN(nSrc, nFire);
alfa     = NaN(nSrc, nFire);
beta     = NaN(nSrc, nFire);
betaErr  = NaN(nSrc, nFire);
pBeta    = NaN(nSrc, nFire);
R2       = NaN(nSrc, nFire);
for i = 1:nSrc
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{i}));
    Size = [data.size]';
    Fire = {data.FireType}';
    for j = 1:nFire
        sz = Size(strcmp(Fire, fireType{j}));
        if isempty(sz), continue; end
        [counts, binEdges] = histcounts(sz, edges);
        binWidth = diff(binEdges);
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
    
        x = log10(binCenters(probDensity > 0)); 
        y = log10(probDensity(probDensity > 0));        
        writetable(table(x', y', 'VariableNames', {'binCenter', 'probDensity'}), ...
            sprintf('dataFig/curve/%s-%s.csv', dataSrc{i}, fireType{j}));

        mdl = fitlm(x, y);
        x_fit = linspace(min(x), max(x), 100)';
        y_fit = predict(mdl, x_fit);
        writetable(table(x_fit, y_fit, 'VariableNames', {'x_fit', 'y_fit'}), ...
            sprintf('dataFig/curve/%s-%s-fit.csv', dataSrc{i}, fireType{j}));

        count(i, j)   = numel(sz);
        area(i, j)    = sum(sz);
        alfa(i, j)    = mdl.Coefficients.Estimate(1);
        beta(i, j)    = mdl.Coefficients.Estimate(2);
        betaErr(i, j) = mdl.Coefficients.SE(2);
        pBeta(i, j)   = mdl.Coefficients.pValue(2);
        R2(i, j)      = mdl.Rsquared.Ordinary;            
    end
end

for i = 1:nSrc
    countVec = squeeze(count(i, :))';
    areaVec = squeeze(area(i, :))';
    alfaVec = squeeze(alfa(i, :))';
    betaVec = squeeze(beta(i, :))';
    betaErrVec = squeeze(betaErr(i, :))';
    pBetaVec = squeeze(pBeta(i, :))';
    R2Vec = squeeze(R2(i, :))';

    writetable(table(countVec, areaVec, alfaVec, betaVec, betaErrVec, pBetaVec, R2Vec, ...
        'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2'}, ...
        'RowNames', fireType), sprintf('dataFig/beta/2Fires-%s.csv',dataSrc{i}), 'WriteRowNames', true);
end
%% beta for 4-datasets, 4-decades, 2-fireTypes
clear; clc;

dataSrc  = {'CalFire', 'MTBS', 'Atlas', 'FIRED'};   nSrc  = numel(dataSrc);
decades  = {'1990s', '2000s', '2010s', '2020s'};    nDec  = numel(decades);
fireType = {'Urban-edge','Wildland'};               nFire = numel(fireType);
count    = NaN(nSrc, nDec, nFire);
area     = NaN(nSrc, nDec, nFire);
alfa     = NaN(nSrc, nDec, nFire);
beta     = NaN(nSrc, nDec, nFire);
betaErr  = NaN(nSrc, nDec, nFire);
pBeta    = NaN(nSrc, nDec, nFire);
R2       = NaN(nSrc, nDec, nFire);


edges = 10 .^ (-1:0.05:5); 
for i = 1:nSrc
    data = shaperead(sprintf('dataPrc/firePrmt/%s.shp', dataSrc{i}));
    Size = [data.size]';
    Fire = {data.FireType}';
    Dec  = {data.decade}';
    for d = 1:nDec
        for j = 1:nFire
            sz = Size(strcmp(Fire, fireType{j}) & strcmp(Dec, decades{d}));
            if isempty(sz), continue; end
            [counts, binEdges] = histcounts(sz, edges);
            binWidth = diff(binEdges);
            probDensity = counts ./ (sum(counts) * binWidth); 
            binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
        
            x = log10(binCenters(probDensity > 0)); 
            y = log10(probDensity(probDensity > 0));
            writetable(table(x', y', 'VariableNames', {'binCenter', 'probDensity'}), ...
                sprintf('dataFig/curve/%s-%s-%s.csv', dataSrc{i}, fireType{j}, decades{d}));

            mdl = fitlm(x, y);
            x_fit = linspace(min(x), max(x), 100)';
            y_fit = predict(mdl, x_fit);
            writetable(table(x_fit, y_fit, 'VariableNames', {'x_fit', 'y_fit'}), ...
                sprintf('dataFig/curve/%s-%s-%s-fit.csv', dataSrc{i}, fireType{j}, decades{d}));

            count(i, d, j)   = numel(sz);
            area(i, d, j)  = sum(sz);
            alfa(i, d, j)   = mdl.Coefficients.Estimate(1);
            beta(i, d, j)   = mdl.Coefficients.Estimate(2);
            betaErr(i, d, j) = mdl.Coefficients.SE(2);
            pBeta(i, d, j)   = mdl.Coefficients.pValue(2);
            R2(i, d, j)      = mdl.Rsquared.Ordinary;            
        end
    end
end

for i = 1:nSrc
    for j = 1:nFire
        countVec = squeeze(count(i, :, j))';
        areaVec = squeeze(area(i, :, j))';
        alfaVec = squeeze(alfa(i, :, j))';
        betaVec = squeeze(beta(i, :, j))'; 
        betaErrVec = squeeze(betaErr(i, :, j))';
        pBetaVec = squeeze(pBeta(i, :, j))';
        R2Vec = squeeze(R2(i, :, j))';

        writetable(table(countVec, areaVec, alfaVec, betaVec, betaErrVec, pBetaVec, R2Vec, ...
            'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2'}, ...
            'RowNames', decades), ...
            sprintf('dataFig/beta/4Decs-%s-%s.csv', dataSrc{i}, fireType{j}), "WriteRowNames", true);
    end
end
%% beta for CalFire, 2-fireTypes, 2-ignition
clear; clc;

fireType = {'Urban-edge','Wildland'}; nFire = numel(fireType);
igType   = {'Human', 'Natural'};      nIg   = numel(igType);

data = shaperead('dataPrc/firePrmt/CalFire.shp');
Size = [data.size]';
Fire = {data.FireType}';
Ign  = {data.Ignition}';

count    = NaN(nFire ,nIg);
area     = NaN(nFire ,nIg);
alfa     = NaN(nFire ,nIg);
beta     = NaN(nFire ,nIg);
betaErr  = NaN(nFire ,nIg);
pBeta    = NaN(nFire ,nIg);
R2       = NaN(nFire ,nIg);
edges = 10 .^ (0:0.05:5); 
for i = 1:nFire
    for j = 1:nIg
        sz = Size(strcmp(Fire, fireType{i}) & strcmp(Ign, igType{j}));
        if isempty(sz), continue; end
        [counts, binEdges] = histcounts(sz, edges);
        binWidth = diff(binEdges);
        probDensity = counts ./ (sum(counts) * binWidth); 
        binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
    
        x = log10(binCenters(probDensity > 0)); 
        y = log10(probDensity(probDensity > 0));
        mdl = fitlm(x, y);

        count(i,j)   = numel(sz);
        area(i,j)    = sum(sz);
        alfa(i,j)    = mdl.Coefficients.Estimate(1);
        beta(i,j)    = mdl.Coefficients.Estimate(2);
        betaErr(i,j) = mdl.Coefficients.SE(2);
        pBeta(i,j)   = mdl.Coefficients.pValue(2);
        R2(i,j)      = mdl.Rsquared.Ordinary;            
    end
end

for i = 1:nFire
    countVec = count(i,:)';  
    areaVec = area(i,:)';    
    alfaVec = alfa(i,:)';    
    betaVec = beta(i,:)';    
    betaErrVec = betaErr(i,:)'; 
    pBetaVec = pBeta(i,:)';   
    R2Vec = R2(i,:)';         

    writetable(table(countVec, areaVec, alfaVec, betaVec, betaErrVec, pBetaVec, R2Vec, ...
        'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2'}, ...
        'RowNames', igType), sprintf('dataFig/beta/2Igns-%s.csv', fireType{i}), "WriteRowNames", true);
end
%% beta for CalFire, four cut-off sizes
edgeCut  = {'<10km2Fire', '<100km2Fire', '<1000km2Fire', 'allFire'}; 
nEdge    = numel(edgeCut);
edges = cell(nEdge,1); 
for k = 1:nEdge, edges{k} = 10.^(0:0.05:k); end

fireType = {'Urban-edge','Wildland'};               nFire = numel(fireType);
decades  = {'1990s', '2000s', '2010s', '2020s'};    nDec  = numel(decades);

data = shaperead('dataPrc/firePrmt/CalFire.shp');
Size = [data.size]';
Fire = {data.FireType}';
Dec  = {data.decade}';

count    = NaN(nEdge, nFire, nDec);
area     = NaN(nEdge, nFire, nDec);
alfa     = NaN(nEdge, nFire, nDec);
beta     = NaN(nEdge, nFire, nDec);
betaErr  = NaN(nEdge, nFire, nDec);
pBeta    = NaN(nEdge, nFire, nDec);
R2       = NaN(nEdge, nFire, nDec);
for i = 1:nEdge
    edge = edges{i};
    for j = 1:nFire
        for d = 1:nDec
            sz = Size(strcmp(Fire, fireType{j}) & strcmp(Dec, decades{d}));   
            if numel(sz) < 10, continue; end
            [counts, binEdges] = histcounts(sz, edge); 
            binWidth = diff(binEdges);
            probDensity = counts ./ (sum(counts) * binWidth); 
            binCenters = sqrt(binEdges(1:end-1) .* binEdges(2:end));
              
            x = log10(binCenters(probDensity > 0)); 
            y = log10(probDensity(probDensity > 0));
            mdl = fitlm(x, y); 
    
            count(i,j,d)   = numel(sz);
            area(i,j,d)    = sum(sz);
            alfa(i,j,d)    = mdl.Coefficients.Estimate(1);
            beta(i,j,d)    = mdl.Coefficients.Estimate(2);
            betaErr(i,j,d) = mdl.Coefficients.SE(2);
            pBeta(i,j,d)   = mdl.Coefficients.pValue(2);
            R2(i,j,d)      = mdl.Rsquared.Ordinary;       
        end    
    end
end

for j = 1:nFire
    for d = 1:nDec
        countVec = squeeze(count(:,j,d));  
        areaVec = squeeze(area(:,j,d));      
        alfaVec = squeeze(alfa(:,j,d));     
        betaVec = squeeze(beta(:,j,d));      
        betaErrVec = squeeze(betaErr(:,j,d));  
        pBetaVec = squeeze(pBeta(:,j,d));  
        R2Vec = squeeze(R2(:,j,d));           
    
        writetable(table(countVec, areaVec, alfaVec, betaVec, betaErrVec, pBetaVec, R2Vec, ...
            'VariableNames', {'count', 'area', 'alfa', 'beta', 'betaErr', 'pBeta', 'R2'}, ...
            'RowNames', edgeCut), sprintf('dataFig/beta/4cutSizes-%s-%s.csv',fireType{j}, decades{d}), 'WriteRowNames', true);
    end
end

%% climate, CalFire
varNames = {'tmean', 'vpdmax','tmin', 'tmax'}; nvar = numel(varNames);
fireType = {'Urban-edge','Wildland'};          nFire = numel(fireType);

data = shaperead('dataPrc/firePrmt/CalFire.shp');
Fire = {data.FireType}';

data_u = data(strcmp(Fire, fireType{1}));
data_w = data(strcmp(Fire, fireType{2}));

Var_u = NaN(numel(data_u), nvar);
Var_w = NaN(numel(data_w), nvar);
for i =1:nvar
    Var_u(:,i) = [data_u.(varNames{i})]';
    Var_w(:,i) = [data_w.(varNames{i})]';
end
writetable(array2table(Var_u, 'VariableNames', string(varNames)),'dataFig/climate/Urban-edge.csv');
writetable(array2table(Var_w, 'VariableNames', string(varNames)),'dataFig/climate/Wildland.csv');
%% NDVI, CalFire
ndviAll = [data.ndviM]';
for i = 1:nFire
    for d = 2:nDec
        ndvi = ndviAll(strcmp(Fire, fireType{i}) & strcmp(Dec, decades{d}) & ~isnan(ndviAll));
        writetable(table(ndvi), sprintf('dataFig/ndvi/%s-%s.csv', fireType{i}, decades{d}));
    end
end
%% year: fire number & burned area, CalFire
data = shaperead('dataPrc/firePrmt/CalFire_fullSize.shp');
Size = [data.size]';
Fire = {data.FireType}';
Ign  = {data.Ignition}';
year = [data.year]';

count    = NaN(nIg, nFire, nyr);
area     = NaN(nIg, nFire, nyr);
for i = 1:nIg
    for j = 1:nFire
        for yr = 1:nyr
            sz = Size(strcmp(Ign, igType{i}) & strcmp(Fire, fireType{j}) & year == years(yr));        
            count(i,j,yr)   = numel(sz);
            area(i,j,yr)    = sum(sz);
        end
    end
end

for j = 1:nFire
    countVec = squeeze(count(:,j,:))'; 
    writetable(array2table(countVec, 'VariableNames', {'Human', 'Natural'}, 'RowNames',string(years)), ...
        sprintf('dataFig/ba/count-%s.csv', fireType{j}), "WriteRowNames",true);

    areaVec = squeeze(area(:,j,:))'; 
    writetable(array2table(areaVec, 'VariableNames', {'Human', 'Natural'}, 'RowNames',string(years)), ...
        sprintf('dataFig/CtBA/area-%s.csv', fireType{j}), "WriteRowNames",true);
end

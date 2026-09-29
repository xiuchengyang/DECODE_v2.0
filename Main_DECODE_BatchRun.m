clear; clc;
close all;

%% Download TMD functions and EOT20 tide model
pathMain = fileparts(mfilename('fullpath'));
addpath(pathMain);
pathTMD = DownloadTMDFunctions(pathMain);
addpath(genpath(pathTMD));
tideModel = DownloadEOT20Model(pathMain);

%% Paths
pathObs = fullfile(pathMain, 'Examples', 'SatelliteObs');
pathOutput = fullfile(pathMain, 'Examples', 'Output');
pathSample = fullfile(pathMain, 'Examples', 'DECODETest.csv');
pathTideLookup = fullfile(pathMain, 'Tide_Lookup.csv');

if ~exist(pathOutput, 'dir'), mkdir(pathOutput); end

%% Settings
sampleIDs = 1:20;
detectBands = {'B2','B3','B4','B8','B11','B12', ...
    'NDVI','EVI','AWEI','MNDWI'};
plotBands = {'B8','B11','NDVI','MNDWI'};
plotLabels = {'NIR','SWIR1','NDVI','MNDWI'};

scale = 10000;
comp = 12;
T_cg = 0.99;
Tmax_cg = 1 - 1e-5;
conse = 6;
startDate = datetime(2017,1,1);
endDate = datetime(2025,12,31);

%% Inputs
sampleInfo = readtable(pathSample, 'VariableNamingRule', 'preserve');
tideLookup = readtable(pathTideLookup, 'VariableNamingRule', 'preserve');

sampleInfo.lat01 = round(sampleInfo.lat * 10) / 10;
sampleInfo.lon01 = round(sampleInfo.lon * 10) / 10;

%% Run samples
for iSample = 1:numel(sampleIDs)

    plotID = sampleIDs(iSample);
    pathObservation = fullfile(pathObs, sprintf('Obs_%03d.csv', plotID));
    pathMat = fullfile(pathOutput, sprintf('DECODE_P%03d.mat', plotID));
    pathFigure = fullfile(pathOutput, ...
        sprintf('P%03d_Landsat_DECODE_TimeSeries.png', plotID));

    idxSample = find(sampleInfo.id == plotID, 1, 'first');
    if isempty(idxSample)
        error('plotID %03d is missing from %s.', plotID, pathSample);
    end

    sampleLat = sampleInfo.lat(idxSample);
    sampleLon = sampleInfo.lon(idxSample);
    sampleLat01 = sampleInfo.lat01(idxSample);
    sampleLon01 = sampleInfo.lon01(idxSample);

    idxTide = find(tideLookup.lat == sampleLat01 & ...
        tideLookup.lon == sampleLon01, 1, 'first');

    if isempty(idxTide)
        tideLat = NaN;
        tideLon = NaN;
        tideFlag = "not used";
    else
        tideLat = tideLookup.use_lat(idxTide);
        tideLon = tideLookup.use_lon(idxTide);
        tideFlag = lower(strtrim(string(tideLookup.flag(idxTide))));
        if ~ismember(tideFlag, ["valid", "nearby"])
            error('Unknown TMD flag for plotID %03d: %s', plotID, tideFlag);
        end
    end

    %% Run DECODE and save MAT
    fprintf('Run DECODE: %03d\n', plotID);

    obs = readtable(pathObservation, 'VariableNamingRule', 'preserve');
    obs = prepareLandsatObservations(obs, scale);

    if height(obs) < conse * 2
        error('Too few observations for plotID %03d.', plotID);
    end

    obs = compositeClearObservations(obs, comp);

    if height(obs) < conse * 2
        error('Too few observations after compositing for plotID %03d.', plotID);
    end

    sdate = obs.sdate;
    line_t = obs{:, detectBands};
    qa = zeros(height(obs), 1);

    if tideFlag == "not used"
        tide = [];
        num_c = 8;
    else
        tide = tmd_predict( ...
            tideModel, tideLat, tideLon, obs.obsDateTime) * 1000;
        tide = double(tide(:));
        tide(isnan(tide)) = 0;
        num_c = 9;
    end

    rec_cg = DECODEPixel( ...
        sdate, line_t, qa, plotID, ...
        T_cg, Tmax_cg, conse, num_c, ...
        'tide', tide, ...
        'numBandGreen', 2, ...
        'numBandSWIR1', 5, ...
        'bandNames', detectBands);

    if isempty(rec_cg)
        error('DECODE returned an empty result for plotID %03d.', plotID);
    end

    save(pathMat, 'rec_cg', 'obs', 'plotID', ...
        'tideLat', 'tideLon', 'tideFlag', ...
        'detectBands', 'scale', '-v7.3');

    %% Plot time series
    if ~all(ismember(plotBands, detectBands))
        error('The DECODE result for plotID %03d does not contain all plot bands.', ...
            plotID);
    end

    fig = figure('Visible', 'on', 'Color', 'w', ...
        'Position', [100 100 1200 950]);
    layout = tiledlayout(fig, 4, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    plotAxes = gobjects(numel(plotBands), 1);

    for iBand = 1:numel(plotBands)
        bandName = plotBands{iBand};
        idxBand = find(strcmp(detectBands, bandName), 1);

        ax = nexttile(layout);
        plotAxes(iBand) = ax;
        hold(ax, 'on');

        hObs = scatter(ax, obs.obsDateTime, obs.(bandName) / scale, ...
            22, [0.25 0.25 0.25], 'filled');
        hFit = gobjects(0);
        hBreak = gobjects(0);

        for iSeg = 1:numel(rec_cg)
            tStart = double(rec_cg(iSeg).t_start(:));
            tEnd = double(rec_cg(iSeg).t_end(:));
            tStart = tStart(isfinite(tStart));
            tEnd = tEnd(isfinite(tEnd));

            if isempty(tStart) || isempty(tEnd)
                error('Invalid DECODE dates for plotID %03d, segment %d.', ...
                    plotID, iSeg);
            end

            xStart = max(min(tStart), datenum(startDate));
            xEnd = min(max(tEnd), datenum(endDate));
            xPlot = (ceil(xStart):floor(xEnd))';

            if ~isempty(xPlot)
                datePlot = datetime(xPlot, 'ConvertFrom', 'datenum');

                if tideFlag == "not used"
                    pred = autoTSPred(xPlot, ...
                        rec_cg(iSeg).coefs(:,idxBand));
                else
                    tideFit = tmd_predict( ...
                        tideModel, tideLat, tideLon, datePlot) * 1000;
                    tideFit = double(tideFit(:));
                    tideFit(isnan(tideFit)) = 0;
                    pred = autoTSPredWL(xPlot, ...
                        rec_cg(iSeg).coefs(:,idxBand), tideFit);
                end

                h = plot(ax, datePlot, pred / scale, ...
                    '-', 'Color', [0.10 0.45 0.80], 'LineWidth', 2);

                if isempty(hFit)
                    hFit = h;
                else
                    h.HandleVisibility = 'off';
                end
            end

            hasChange = any(double(rec_cg(iSeg).change_prob(:)) == 1);
            breakNumbers = double(rec_cg(iSeg).t_break(:));
            breakNumbers = breakNumbers(isfinite(breakNumbers));

            if hasChange && ~isempty(breakNumbers)
                data = rec_cg(iSeg).data;

                for iBreak = 1:numel(breakNumbers)
                    breakDate = datetime(breakNumbers(iBreak), ...
                        'ConvertFrom', 'datenum');
                    [~, idxNear] = min(abs( ...
                        double(data.sdate) - breakNumbers(iBreak)));
                    breakValue = double(data.(bandName)(idxNear)) / scale;

                    h = scatter(ax, breakDate, breakValue, 90, ...
                        'MarkerEdgeColor', [0.90 0.30 0.10], ...
                        'MarkerFaceColor', 'none', 'LineWidth', 1.6);

                    if isempty(hBreak)
                        hBreak = h;
                    else
                        h.HandleVisibility = 'off';
                    end
                end
            end
        end

        if isempty(hFit)
            error('No fitted %s curve for plotID %03d.', ...
                bandName, plotID);
        end

        xlim(ax, [startDate, datetime(2026,1,1)]);
        xticks(ax, datetime(2017:2026,1,1));
        xtickformat(ax, 'yyyy');
        ylabel(ax, plotLabels{iBand});
        grid(ax, 'on');
        box(ax, 'on');
        set(ax, 'FontName', 'Arial', 'FontSize', 13);

        if ismember(bandName, {'B8', 'B11'})
            yLimits = ylim(ax);
            ylim(ax, [0, yLimits(2)]);
        elseif ismember(bandName, {'NDVI', 'MNDWI'})
            ylim(ax, [-1 1]);
        end

        if iBand < numel(plotBands)
            ax.XTickLabel = [];
        else
            xlabel(ax, 'Year');
        end

        if iBand == 1
            legendHandles = [hObs, hFit];
            legendNames = {'Landsat observations', ...
                'DECODE fitted time series'};
            if ~isempty(hBreak)
                legendHandles = [legendHandles, hBreak];
                legendNames{end+1} = 'Change point';
            end
            legend(ax, legendHandles, legendNames, ...
                'Location', 'best');
        end
    end

    linkaxes(plotAxes, 'x');
    title(layout, sprintf( ...
        'Pixel %03d | Lat/Lon: %.6f, %.6f | TMD: %s', ...
        plotID, sampleLat, sampleLon, tideFlag), ...
        'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 15);

    exportgraphics(fig, pathFigure, 'Resolution', 300);

    fprintf('Saved: %s\n', pathFigure);
end

fprintf('Done.\n');

%% Functions
function T = prepareLandsatObservations(T, outputScale)

    if ismember('fmask', T.Properties.VariableNames)
        T = T(ismember(T.fmask, [0 1]), :);
    end

    dateText = string(T.date) + " " + string(T.time);
    T.obsDateTime = datetime(dateText, ...
        'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    T.date = dateshift(T.obsDateTime, 'start', 'day');
    T.sdate = datenum(T.obsDateTime);

    B2 = double(T.blue);
    B3 = double(T.green);
    B4 = double(T.red);
    B8 = double(T.nir);
    B11 = double(T.swir1);
    B12 = double(T.swir2);

    T.B2 = B2 * outputScale;
    T.B3 = B3 * outputScale;
    T.B4 = B4 * outputScale;
    T.B8 = B8 * outputScale;
    T.B11 = B11 * outputScale;
    T.B12 = B12 * outputScale;

    T.NDVI = safeDivide(B8 - B4, B8 + B4) * outputScale;
    T.EVI = 2.5 * safeDivide( ...
        B8 - B4, B8 + 6 * B4 - 7.5 * B2 + 1) * outputScale;
    T.AWEI = (B2 + 2.5 * B3 - 1.5 * (B8 + B11) ...
        - 0.25 * B12) * outputScale;
    T.MNDWI = safeDivide(B3 - B11, B3 + B11) * outputScale;

    keep = isfinite(T.sdate) & isfinite(T.NDVI) & ...
        isfinite(T.EVI) & isfinite(T.AWEI) & isfinite(T.MNDWI);
    T = sortrows(T(keep, :), 'obsDateTime');
end

function y = safeDivide(a, b)
    y = zeros(size(a));
    idx = isfinite(a) & isfinite(b) & b ~= 0;
    y(idx) = a(idx) ./ b(idx);
end

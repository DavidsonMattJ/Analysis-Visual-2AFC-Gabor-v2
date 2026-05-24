function plot_GaitresultsBinned_2speed_vAUD(dataIN, GFX_FourierNull, cfg)
% plot_GaitresultsBinned_2speed_vAUD
%
% Plots gait-phase-resolved SDT time-series (within-gait) for the auditory
% detection task, at either group (GFX) or participant (PFX) level.
%
% Each subplot shows a chosen DV (e.g. d', HR, Accuracy) as a function of
% position within the gait cycle (0–100%), with optional overlays for the
% best-fit Fourier curve, Fourier R²(Hz) profile, and head-height trace.
%
% Arguments
% ---------
%   dataIN          [nSubs × nSpeeds × nFoot] struct — binned SDT data
%                   (GFX_TargPosData or GFX_RespPosData; j4b output)
%   GFX_FourierNull [1 × nSpeeds] struct — group-level R²(Hz) profiles (j6 output)
%                   Pass [] when unavailable; disables plotFITforced / plotShuff.
%   PFX_FourierNull [nSubs × nSpeeds] struct — participant-level R²(Hz) profiles
%                   (j9 output). Reserved for future use; pass [] for now.
%   cfg             Configuration struct — see fields below.
%
% cfg fields
% ----------
%   plotlevel     'GFX' (group mean ± SEM) | 'PFX' (one figure per participant)
%   usebin        1 = binned data (recommended); 0 = raw 100-point gait% series
%   type          'Target' | 'Response'  (onset type, used in labels and field names)
%   DV            'Accuracy' | 'RT' | 'dprime' | 'crit' | 'counts' | 'FA' | 'HR'
%   binDV         'Calc' (pooled-bin SDT recalculation) | 'Av' (within-bin mean)
%   nGaitstoplot  1 = single step (gc); 2 = stride (doubgc); [1 2] = both
%   ispeed        speed indices to plot: 1=slow, 2=normal, 3=combined
%   plotCOL       data line/marker colour
%   fitCOL        best-fit line colour (also tints the DV ylabel)
%   yyaxis        'left' | 'right'  (which y-axis to use for the DV)
%   normON        1 = normalise per participant before averaging; 0 = raw values
%   normtype      'absolute'|'relative'|'relchange'|'normchange'|'db'
%   plotHead      1 = overlay mean head-height trace on right y-axis
%   plotFIT       1 = overlay unbounded fourier1 best fit on the data panel
%   plotFITforced 1 = show Fourier R²(Hz) frequency-profile panel (needs GFX_FourierNull)
%   plotShuff     1 = overlay null 95th-percentile CV on R²(Hz) panel
%   plotGaitGrid  1 = 3×3 grid (all foot × speed); 0 = single line per speed
%   plotnewLayout 1 = focused 3-column layout [LR | RL | LRL], one row per speed
%   omitPpants    exclusion indices (applied upstream in j7; leave [] here)
%   pidx1         bin-edge vector for single-step (gc) gait cycle
%   pidx2         bin-edge vector for stride (doubgc) gait cycle
%   HeadData      GFX_headY array [nSubs × nSpeeds]
%   subjIDs       cell array of participant ID strings
%   figdir        root output directory for PNG saves
%   datadir       base data directory (passed through; not used internally)

% -------------------------------------------------------------------------
%  Unpack frequently used cfg fields
% -------------------------------------------------------------------------
GFX_headY = cfg.HeadData;
nsubs      = length(cfg.subjIDs);
pcntindexes = {cfg.pidx1, cfg.pidx2};

% Colour palettes (consistent across all layouts)
usecols   = {[0 .7 0], [.7 0 0], [.7 0 .7]};       % foot: L=green, R=red, comb=purple
speedCols = {'b', [1, 171/255, 64/255], 'k'};        % speed: slow=blue, normal=amber, comb=black

% String labels
ftnames_step   = {'LR',  'RL',  'combined'};         % single-step foot labels
ftnames_stride = {'LRL', 'RLR', 'combined'};         % stride foot labels
speednames     = {'Slow', 'Normal', 'Combined'};
cyclenames     = {'step', 'stride'};
binnames       = {'unbinned', 'binned'};
usegaitnames   = {'gc', 'doubgc'};

% -------------------------------------------------------------------------
%  Build the struct-field suffix for the requested DV × binning combination.
%
%  cfg.binDV ('Calc' or 'Av') selects between:
%    Calc – pooled-bin SDT recalculation (preferred; matches j4b / j6 / j9)
%    Av   – simple within-bin mean (less preferred for SDT metrics)
%
%  Field names mirror j4b output, e.g.:  doubgc_binned_AccCalc
%                                         gc_binned_HRCalc
% -------------------------------------------------------------------------
if cfg.usebin
    switch cfg.DV
        case 'RT',       usefield = ['_binned_rts'  cfg.binDV];
        case 'Accuracy', usefield = ['_binned_Acc'  cfg.binDV];
        case 'crit',     usefield = ['_binned_crit' cfg.binDV];
        otherwise        % dprime, HR, FA, counts — field names are consistent
            usefield = ['_binned_' cfg.DV cfg.binDV];
    end
else
    % Unbinned: raw 100-point gait-% series (no Calc/Av suffix)
    switch cfg.DV
        case 'RT',       usefield = '_rts';
        case 'Accuracy', usefield = '_Acc';
        otherwise,       usefield = ['_' cfg.DV];
    end
end

% -------------------------------------------------------------------------
%  Figure layout: add a rightmost column for the Fourier R²(Hz) panel when
%  plotFITforced = 1; otherwise use 3 columns (one per speed).
%  Only relevant for GFX single-line layout.
% -------------------------------------------------------------------------
if cfg.plotFITforced && cfg.usebin
    nCols = 4;   % [speed1 | speed2 | speed3 | Fourier R²]
else
    nCols = 3;
end

% Open a figure window 
figure(1); clf;
set(gcf, 'Color', 'w', 'Units', 'normalized', 'Position', [0 0 .6 .9]);

% Default gait type used in filenames (overridden inside loops below)
nGaits_toPlot = cfg.nGaitstoplot(end);



% =========================================================================
%%  GFX — group-level figure
% =========================================================================
if strcmp(cfg.plotlevel, 'GFX')

    psubj = 'GFX';

    % -----------------------------------------------------------------
    %%  GFX standard layout (single-line per speed, or 3×3 grid)
    % -----------------------------------------------------------------
    if ~cfg.plotnewLayout

        for nGaits_toPlot = cfg.nGaitstoplot

            pidx    = pcntindexes{nGaits_toPlot};
            xvec    = computeXvec(pidx, cfg.usebin);
            ftnames = selectFtnames(nGaits_toPlot, ftnames_step, ftnames_stride);

            % ---------------------------------------------------------
            %%  GFX single-line: combined foot, one subplot per speed
            % ---------------------------------------------------------
            if ~cfg.plotGaitGrid

                for iSpeed = cfg.ispeed
                    iLR = 3;   % combined foot only (L/R split in grid mode below)

                    % Collate data across participants
                    ppantData = zeros(nsubs, length(xvec));
                    headData  = [];
                    for isub = 1:nsubs
                        ppantData(isub, :) = dataIN(isub, iSpeed, iLR).( ...
                            [usegaitnames{nGaits_toPlot} usefield]);
                        headData(isub, :)  = GFX_headY(isub, iSpeed).( ...   %#ok<AGROW>
                            [usegaitnames{nGaits_toPlot}]);
                    end

                    % Optional per-participant normalisation
                    if cfg.normON
                        ppantData = applyNorm(ppantData, cfg.normtype);
                         % Normalisation suffix appended to output field names ('norm' or empty).
                        ntype = 'norm';
                    else
                        ntype = [];

                    end

                    % Replace Inf with NaN before group averaging
                    ppantData(isinf(ppantData)) = nan;
                    gM  = mean(ppantData, 1, 'omitnan');
                    stE = CousineauSEM(ppantData);

                    if all(isnan(gM))
                        fprintf('  Skipping iSpeed=%d for %s — all values NaN.\n', iSpeed, cfg.DV);
                        continue
                    end

                    % Data panel
                    ax = subplot(2, nCols, iSpeed+ nCols*(nGaits_toPlot-1));
                    hold(ax, 'on');

                    if strcmp(cfg.yyaxis, 'right')
                        yyaxis(ax, 'right');
                    end

                    plot(ax, xvec, gM, 'o-', 'Color', cfg.plotCOL, 'LineWidth', 2);
                    errorbar(ax, xvec, gM, stE, ...
                        'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2);

                    % Auto-scale to ±1× the data range to capture variance
                    sdrange = max(gM) - min(gM);
                    ylim(ax, [min(gM) - sdrange, max(gM) + sdrange]);
                    box(ax, 'on');

                    % Unbounded fourier1 best fit (a0 + a1·cos(w·x) + b1·sin(w·x))
                    if cfg.plotFIT
                        [f, gof] = fit(xvec', gM', 'fourier1');
                        h = plot(f, xvec, gM);
                        h(2).LineWidth = 4;
                        h(2).Color     = speedCols{iSpeed};
                        Hzapp      = xvec(end) / (2*pi / f.w);
                        legdetails = [sprintf('%.2f', Hzapp) ' Hz_{GC},  R^2 = ' ...
                            sprintf('%.2f', gof.rsquare)];
                        legend(h(2), legdetails, 'FontSize', 15, ...
                            'AutoUpdate', 'off', 'Location', 'best');
                    end

                    % Head-height trace overlay (normalised, on right y-axis)
                    if cfg.plotHead
                        headMean = mean(headData, 1, 'omitnan');
                        headMean = headMean ./ max(headMean);              % normalise 0–1
                        yl       = ylim(ax);
                        headMean = headMean .* (yl(2) + 0.1*diff(yl));    % scale above data
                        headMean = imresize(headMean, [1, length(xvec)]);
                        stEH     = imresize(CousineauSEM(headData), [1, length(xvec)]);
                        yyaxis(ax, 'right');
                        sh = shadedErrorBar(xvec, headMean, stEH, ...
                            {'Color', speedCols{iSpeed}}, 1);
                        sh.mainLine.LineWidth = 2;
                        set(ax, 'YTick', []);
                        yl2 = ylim(ax);
                        ylim(ax, [yl2(1) - 1.5*diff(yl2), yl2(2)]);  % drop floor to show data
                    end

                    title(ax, [psubj '  ' speednames{iSpeed} '  N=' num2str(nsubs)], ...
                        'Interpreter', 'none');
                    xlabel(ax, [cfg.type ' position as % of ' cyclenames{nGaits_toPlot} '-cycle']);
                    ylabel(ax, cfg.DV, 'FontWeight', 'bold', 'Color', cfg.fitCOL);
                    set(ax, 'YColor', 'k');
                    formatGaitAxis(ax, xvec);

                    % -------------------------------------------------
                    %%  Fourier R²(Hz) panel — rightmost column
                    %   All speeds are overlaid on a single panel using
                    %   'hold on' across iSpeed iterations.
                    %   Solid line = observed R²; dashed = null 95th pct.
                    % -------------------------------------------------
                    if cfg.plotFITforced && cfg.usebin && ~isempty(GFX_FourierNull)

                        Hzspace = 0.01:0.2:10;   % must match j6 frequency grid
                        axF = subplot(2, nCols, nCols+ nCols*(nGaits_toPlot-1));
                        hold(axF, 'on');

                        fldObs  = [cfg.type 'Ons_' usegaitnames{nGaits_toPlot} ...
                            usefield '_fitsRsq_Obs' ntype];
                        fitsObs = GFX_FourierNull(iSpeed).(fldObs);

                        plot(axF, Hzspace(1:length(fitsObs)), fitsObs, ...
                            'Color', speedCols{iSpeed}, 'LineWidth', 3);

                        if cfg.plotShuff
                            % Row 3 of ShuffCV = 95th-percentile null R²
                            fldCV  = [cfg.type 'Ons_' usegaitnames{nGaits_toPlot} ...
                                usefield '_fitsRsq_ShuffCV'];
                            fitsCV = GFX_FourierNull(iSpeed).(fldCV);
                            plot(axF, Hzspace(1:length(fitsObs)), fitsCV(3, :), ...
                                ':', 'LineWidth', 2, 'Color', speedCols{iSpeed});
                        end

                        xlabel(axF, 'Frequency (cycles per stride)');
                        ylabel(axF, 'R^2');
                        ylim(axF, [0 1]);
                        xlim(axF, [0 max(Hzspace)/2])

                        title(axF, 'Fourier R^2 profile', 'FontSize', 15);
                        set(axF, 'FontSize', 15);

                    end

                end % iSpeed (single-line)

            % ---------------------------------------------------------
            %%  GFX 3×3 grid: all foot × speed combinations
            % ---------------------------------------------------------
            else

                pc = 1;

                for iSpeed = 1:3
                    for iLR = 1:3

                        ppantData = zeros(nsubs, length(xvec));
                        headData  = [];
                        for isub = 1:nsubs
                            ppantData(isub, :) = dataIN(isub, iSpeed, iLR).( ...
                                [usegaitnames{nGaits_toPlot} usefield]);
                            headData(isub, :)  = GFX_headY(isub, iSpeed).( ...  %#ok<AGROW>
                                [usegaitnames{nGaits_toPlot}]);
                        end

                        if cfg.normON
                            ppantData = applyNorm(ppantData, cfg.normtype);
                        end

                        gM  = mean(ppantData, 1, 'omitnan');
                        stE = CousineauSEM(ppantData);

                        ax = subplot(3, 3, pc);
                        hold(ax, 'on');
                        yyaxis(ax, 'left');

                        plot(ax, xvec, gM, 'o-', 'Color', cfg.plotCOL, 'LineWidth', 2);
                        errorbar(ax, xvec, gM, stE, ...
                            'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2);
                        ylabel(ax, cfg.DV);

                        sdrange = max(gM) - min(gM);
                        ylim(ax, [min(gM) - sdrange, max(gM) + sdrange]);
                        box(ax, 'on');

                        if cfg.plotFIT
                            [f, gof] = fit(xvec', gM', 'fourier1');
                            h = plot(f, xvec, gM);
                            h(2).LineWidth = 4;
                            h(2).Color     = speedCols{iSpeed};
                            Hzapp      = xvec(end) / (2*pi / f.w);
                            legdetails = [sprintf('%.2f', Hzapp) ' Hz_{GC},  R^2 = ' ...
                                sprintf('%.2f', gof.rsquare)];
                            legend(h(2), legdetails, 'FontSize', 15, ...
                                'AutoUpdate', 'off', 'Location', 'best');
                        end

                        if cfg.plotHead
                            yyaxis(ax, 'right');
                            headMean = imresize(mean(headData, 1, 'omitnan'), [1, xvec(end)]);
                            plot(ax, headMean, 'o', 'Color', speedCols{iSpeed}, 'LineWidth', 3);
                            set(ax, 'YTick', []);
                        end

                        title(ax, [psubj '  ' speednames{iSpeed} '  ' ftnames{iLR}], ...
                            'Interpreter', 'none');
                        xlabel(ax, [cfg.type ' onset (% of ' cyclenames{nGaits_toPlot} '-cycle)']);
                        formatGaitAxis(ax, xvec);
                        pc = pc + 1;

                    end % iLR
                end % iSpeed (grid)

            end % plotGaitGrid

        end % nGaits_toPlot

    % -----------------------------------------------------------------
    %%  GFX focused layout: [LR | RL | LRL] columns, one row per speed
    % -----------------------------------------------------------------
    else

        clf
        % Column assignment: LR (step) | RL (step) | LRL (stride, combined)
        tmpGaits = [1, 1, 2];
        tmpiLR   = [1, 2, 2];
        pc = 1;

        for iSpeed = 1:3
            for istep = 1:3

                nGaits_toPlot = tmpGaits(istep);
                iLR           = tmpiLR(istep);
                pidx          = pcntindexes{nGaits_toPlot};
                xvec          = computeXvec(pidx, cfg.usebin);
                ftnames       = selectFtnames(nGaits_toPlot, ftnames_step, ftnames_stride);

                ppantData = zeros(nsubs, length(xvec));
                headData  = [];
                rawCounts = [];
                for isub = 1:nsubs
                    ppantData(isub, :) = dataIN(isub, iSpeed, iLR).( ...
                        [usegaitnames{nGaits_toPlot} usefield]);
                    headData(isub, :)  = GFX_headY(isub, iSpeed).( ...    %#ok<AGROW>
                        [usegaitnames{nGaits_toPlot}]);
                    rawCounts(isub, :) = dataIN(isub, iSpeed, iLR).( ...  %#ok<AGROW>
                        [usegaitnames{nGaits_toPlot} '_counts']);
                end

                if cfg.normON
                    ppantData = applyNorm(ppantData, cfg.normtype);
                end

                gM          = mean(ppantData, 1, 'omitnan');
                stE         = CousineauSEM(ppantData);
                totalCounts = mean(sum(rawCounts, 2), 'omitnan');

                ax = subplot(3, 3, pc);
                hold(ax, 'on');
                yyaxis(ax, 'left');

                bar(ax, xvec, gM, 'FaceColor', usecols{iLR}, 'FaceAlpha', 0.7);
                errorbar(ax, xvec, gM, stE, ...
                    'Color', 'k', 'LineStyle', 'none', 'LineWidth', 2);
                ylabel(ax, cfg.DV);

                yyaxis(ax, 'right');
                headMean = imresize(mean(headData, 1, 'omitnan'), [1, xvec(end)]);
                plot(ax, headMean, 'o', 'Color', speedCols{iSpeed}, 'LineWidth', 3);
                set(ax, 'YTick', []);

                title(ax, {[psubj '  ' speednames{iSpeed} '  ' ftnames{iLR}]; ...
                    ['raw counts: ' num2str(round(totalCounts))]}, 'Interpreter', 'none');
                xlabel(ax, [cfg.type ' onset (% of ' cyclenames{nGaits_toPlot} '-cycle)']);
                formatGaitAxis(ax, xvec);
                pc = pc + 1;

            end % istep
        end % iSpeed

    end % plotnewLayout

    % Save GFX figure — one PNG per call (last gait type in filename)
    outDir = fullfile(cfg.figdir, [cfg.type ' onset ' cfg.DV ' binned']);
    ensureDir(outDir);
    print(fullfile(outDir, [psubj ' ' cfg.type ' onset ' cfg.DV ' ' ...
        binnames{cfg.usebin+1} ' ' usegaitnames{nGaits_toPlot}]), '-dpng');
    shg

end % GFX
end % main function


% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function xvec = computeXvec(pidx, usebin)
% computeXvec  Compute bin-centre x-axis values from a bin-edge vector.
%   usebin = 1: centre of each bin (standard for binned data).
%   usebin = 0: full 1:pidx(end) range (all 100 gait% points).
if usebin
    mdiff = round(mean(diff(pidx)) / 2);
    xvec  = pidx(1:end-1) + mdiff;
else
    xvec = 1:pidx(end);
end
end


function data = applyNorm(data, normtype)
% applyNorm  Per-participant normalisation of a [nSubs × nBins] matrix.
%
%   Each row is normalised relative to that participant's row mean, so
%   the normalised grand mean across bins is approximately zero (or one).
%
%   normtype options:
%     'absolute'  – subtract row mean
%     'relative'  – divide by row mean, minus 1  (percent change from mean)
%     'relchange' – (x – mean) / mean             (classic relative change)
%     'normchange'– (x – mean) / (x + mean)       (bounded symmetric change)
%     'db'        – 10·log10(x / mean)             (decibel change)
pM       = mean(data, 2, 'omitnan');
meanVals = repmat(pM, 1, size(data, 2));
switch normtype
    case 'absolute',   data = data - meanVals;
    case 'relative',   data = data ./ meanVals - 1;
    case 'relchange',  data = (data - meanVals) ./ meanVals;
    case 'normchange', data = (data - meanVals) ./ (data + meanVals);
    case 'db',         data = 10 * log10(data ./ meanVals);
    otherwise
        warning('applyNorm: unknown normtype ''%s'' — returning raw data.', normtype);
end
end


function formatGaitAxis(ax, xvec)
% formatGaitAxis  Apply standard gait-cycle axis tick formatting.
%   Labels start (0%), midpoint (50%), and end (100%) of the gait cycle.
midp = xvec(ceil(length(xvec) / 2));
set(ax, 'FontSize', 15, ...
    'XTick',      [xvec(1), midp, xvec(end)], ...
    'XTickLabel', {'0%', '50%', '100%'});
end


function ft = selectFtnames(nGaits_toPlot, step_names, stride_names)
% selectFtnames  Return the appropriate foot-condition label cell array.
%   nGaits_toPlot = 1 → single-step labels {LR, RL, combined}
%   nGaits_toPlot = 2 → stride labels      {LRL, RLR, combined}
if nGaits_toPlot == 1
    ft = step_names;
else
    ft = stride_names;
end
end


function ensureDir(dirPath)
% ensureDir  Create dirPath (and any missing parents) if it does not exist.
if ~isfolder(dirPath)
    mkdir(dirPath);
end
end

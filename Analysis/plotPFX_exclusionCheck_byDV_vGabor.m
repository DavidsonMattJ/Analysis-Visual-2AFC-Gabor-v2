% plotPFX_exclusionCheck_byDV.m
%
% Variant of plotPFX_exclusionCheck_targetCounts with two extensions:
%   1. Summary figure shows all three walk speeds as separate columns
%      (Slow | Normal | Combined), rather than combined-only.
%   2. The dependent variable (DV) is configurable at the top of the script,
%      so the same figure logic works for target counts, accuracy, FA, etc.
%
% See plotPFX_exclusionCheck_targetCounts for full documentation of the
% underlying approach (uniformity check, Fourier R² profiles, etc.).
%
% DV CONFIGURATION (set in the block below):
%   typeDV– identifier used in output filenames / figure titles
%   fieldName_step  – field in GFX_TargPosData for the step (gc) binned DV
%   fieldName_stride– field in GFX_TargPosData for the stride (doubgc) binned DV
%   dvYLim          – y-axis limits for bar panels ([] = auto)
%
% Layout of per-participant figures (3 rows × 4 cols):
%   Rows    : walk speed         (Slow | Normal | Combined)
%   Columns : [step DV | step R²(Hz) | stride DV | stride R²(Hz)]
%
% Layout of summary figure (2 rows × 3 cols):
%   Rows    : step R² at FOI | stride R² at FOI
%   Columns : Slow | Normal | Combined
%
% Prerequisites:
%   j4b_gaitPhaseResolved_vAUD   – GFX_Data_inGaits[_new].mat
%   j9 full-perm                 – GFX_Data_inGaits_FourierFits[_new].mat
%     (required: provides _fitsRsq_Obs for R² panels and _fitsRsq_ShuffCV
%      for null CV overlays and significance asterisks)
%
% Saves PNGs to: figdir/Participant exclusion check - <dvName>/

% -------------------------------------------------------------------------
%  Directories (guard: skip if homedir already on path)
% -------------------------------------------------------------------------
if ~contains(path, homedir)
    set_myOnlineDirectories_vGaborv2;
end

% =========================================================================
%  DV CONFIGURATION  — edit this block to switch dependent variable
% =========================================================================
testDVs  = {'Accuracy', 'RT', 'dprime', 'crit', 'HR', 'FA', 'counts', 'counts'};

testtype=7; % .

% Example:
% dvName          = 'Target counts';tr               % used in filenames / titles
% fieldName_step  = 'gc_binned_countsCalc';       % field in GFX_TargPosData (step)
% fieldName_stride= 'doubgc_binned_countsCalc';   % field in GFX_TargPosData (stride)
dvYLim          = [];                           % [] = auto; or e.g. [0 1]

dvName = ['Target ' testDVs{testtype}];


% =========================================================================
%  JOB TOGGLES
% =========================================================================
runParticipantPlots = 1;   % 1 = generate + save one figure per participant
runSummaryPlot      = 0;   % 1 = generate + save group ranking figure

% plotPermCV: overlay null 95th-percentile CV band on per-participant R² panels.
% Requires PFX_FourierNull (j9 full-perm output).
plotPermCV = 1;
overlayFIT =1 ; % overlay unbounded Fourier fit to DV.
% -------------------------------------------------------------------------
%  File-name suffix
% -------------------------------------------------------------------------
if ~exist('usenewTrialStride', 'var')
    usenewTrialStride = 1;
end
if usenewTrialStride
    appendp = '_new';
else
    appendp = '';
end

%% Load data ------------------------------------------------------------------
cd([procdatadir filesep 'GFX'])
load(['GFX_Data_inGaits' appendp], ...
    'GFX_TargPosData', 'pidx1', 'pidx2', 'subjIDs');

% PFX_FourierNull provides both observed R² profiles (_fitsRsq_Obs) and
% null CV bands (_fitsRsq_ShuffCV). Always load; R² panels are skipped if
% unavailable (bar-chart panels still render from GFX_Data_inGaits alone).
hasPFXNull = false;
try
    load(['GFX_Data_inGaits_FourierFits' appendp], 'PFX_FourierNull', 'Hzspace');
    hasPFXNull = true;
catch
    fprintf(['Note: PFX_FourierNull not found in ' ...
        'GFX_Data_inGaits_FourierFits%s.mat\n' ...
        '      Run j9 with testFull_perm=1, testtype=7, nGaits_toPlot=1:2.\n' ...
        '      R² panels, CV overlays, and significance asterisks will be skipped.\n'], appendp);
    plotPermCV = 0;
end

%% Derived x-axis vectors and Nyquist limits ----------------------------------
mdiff1  = round(mean(diff(pidx1)) / 2);
xvec1   = pidx1(1:end-1) + mdiff1;
nyqlim1 = floor(length(xvec1) / 2);
hzlim1 = floor(length(Hzspace)/2);

mdiff2  = round(mean(diff(pidx2)) / 2);
xvec2   = pidx2(1:end-1) + mdiff2;
nyqlim2 = floor(length(xvec2) / 2);
hzlim2 = floor(length(Hzspace)/2);

%% Frequency axis — must match the Hzspace used in j6/j9 ---------------------
% Only the vector is needed here (no fitting); pre-computed R² profiles are
% loaded from PFX_FourierNull (_fitsRsq_Obs / _fitsRsq_ShuffCV).
Hzspace = 0.01:0.2:10;

%% Plot configuration ---------------------------------------------------------
fontsize    = 15;
% stepCol     = [0.00, 0.55, 0.27];
% strideCol   = [0.10, 0.30, 0.80];
speedCols = {'b', [1, 171/255, 64/255], 'k'};
speedLabels = {'Slow', 'Normal', 'Combined'};
iLR         = 3;
ncols       = 4;
nSubs       = size(GFX_TargPosData, 1);

inspectDir = fullfile(figdir, ['Participant exclusion check - ' dvName]);
if ~isfolder(inspectDir)
    mkdir(inspectDir);
end

if runParticipantPlots
    set(gcf, 'Color', 'w', 'Units', 'normalized', 'Position', [0 0 1 1]);
end

%% Frequency of interest ------------------------------------------------------
freqOfInterest_step   = 1;
freqOfInterest_stride = 2;

[~, foiIdx_step]   = min(abs(Hzspace - freqOfInterest_step));
[~, foiIdx_stride] = min(abs(Hzspace - freqOfInterest_stride));
% will also have option for the max (below).

fprintf('Summary FOI — step: %.2f cps (idx %d) | stride: %.2f cps (idx %d)\n', ...
    Hzspace(foiIdx_step), foiIdx_step, Hzspace(foiIdx_stride), foiIdx_stride);

%% Storage for summary figures -------------------------------------------------
% [nSubs × 3 speeds] for both observed R² and null 95th-pct CV
summaryRsq_step      = nan(nSubs, 3); % storing at FOI.
summaryRsq_stride    = nan(nSubs, 3);
summaryNullCV_step   = nan(nSubs, 3);
summaryNullCV_stride = nan(nSubs, 3);

summaryRsqMax_step      = nan(nSubs, 3); % storing at max across CPS
summaryRsqMax_stride    = nan(nSubs, 3);
summaryNullCVMax_step   = nan(nSubs, 3);
summaryNullCVMax_stride = nan(nSubs, 3);
summaryRsqHzidMax_step   = nan(nSubs, 3); % ened extra index for the Hz at max rsq (not FOI).
summaryRsqHzidMax_stride = nan(nSubs, 3);


%% determine fields for look up:
typeDV    = testDVs{testtype};
% Select struct field name for this DV ------------------------------------
% Note: RT uses '_rtsCalc' (the SDT-derived calculated value from j4b),
% consistent with j6. The binned-average field '_rtsAv' is also available
% in the null data (saved by j5) but '_rtsCalc' is preferred for uniformity.
gaitfield = {'gc', 'doubgc'};
binfield  = {'', '_binned'};

if strcmp(typeDV, 'RT')
    fieldName_step = 'gc_binned_rtsCalc';
    fieldName_stride = 'doubgc_binned_rtsCalc';
elseif strcmp(typeDV, 'Accuracy')
    fieldName_step= 'gc_binned_AccCalc';
    fieldName_stride= 'doubgc_binned_AccCalc';
elseif strcmp(typeDV, 'crit')
    fieldName_step= ['gc_binned_critCalc'];
    fieldName_stride= ['doubgc_binned_critCalc'];
else
    % Covers: dprime, HR, FA, counts (field-name pattern is consistent)
    fieldName_step =['gc_binned_' typeDV 'Calc'];
    fieldName_stride =['doubgc_binned_' typeDV 'Calc'];
end


%% Per-participant loop --------------------------------------------------------
for ippant = 1:nSubs

    if runParticipantPlots
        clf
    end

    for ispeed = 1:3

        %% Retrieve DV data for this participant × speed ----------------------
        stepData   = GFX_TargPosData(ippant, ispeed, iLR).(fieldName_step);
        strideData = GFX_TargPosData(ippant, ispeed, iLR).(fieldName_stride);

        %% Look up pre-computed R²(Hz) profiles from PFX_FourierNull ----------
        % Field names are constructed from the DV fieldName variables so this
        % works for any DV that j9 has fitted (e.g. accuracy, FA, counts).
        %   _fitsRsq_Obs    : R² at each Hzspace frequency, observed data
        %   _fitsRsq_ShuffCV: [3 × nFreq] — rows 1/2/3 = .05/.5/.95 quantiles
        rsq_step   = nan(1, length(Hzspace));
        rsq_stride = nan(1, length(Hzspace));
        nullCV_step   = [];
        nullCV_stride = [];

        if hasPFXNull
            try
                rsq_step = PFX_FourierNull(ippant, ispeed). ...
                    (['TargetOns_' fieldName_step '_fitsRsq_Obs']);
            catch
            end
            try
                rsq_stride = PFX_FourierNull(ippant, ispeed). ...
                    (['TargetOns_' fieldName_stride '_fitsRsq_Obs']);
            catch
            end
            try
                nullCV_step = PFX_FourierNull(ippant, ispeed). ...
                    (['TargetOns_' fieldName_step '_fitsRsq_ShuffCV']);
            catch
            end
            try
                nullCV_stride = PFX_FourierNull(ippant, ispeed). ...
                    (['TargetOns_' fieldName_stride '_fitsRsq_ShuffCV']);
            catch
            end
        end

        %% Store observed R² and null 95th-pct CV at FOI ---------------------
        summaryRsq_step  (ippant, ispeed) = rsq_step  (foiIdx_step);
        summaryRsq_stride(ippant, ispeed) = rsq_stride(foiIdx_stride);

        if ~isempty(nullCV_step)
            summaryNullCV_step  (ippant, ispeed) = nullCV_step  (3, foiIdx_step);
        end
        if ~isempty(nullCV_stride)
            summaryNullCV_stride(ippant, ispeed) = nullCV_stride(3, foiIdx_stride);
        end
        % same now for the max values.
        %summaryRsqMax_step      = nan(nSubs, 3); % storing at max across CPS
        % summaryRsqMax_stride    = nan(nSubs, 3);
        % summaryNullCVMax_step   = nan(nSubs, 3);
        % summaryNullCVMax_stride = nan(nSubs, 3);
        % summaryRsqHzidMax_step   = nan(nSubs, 3); % ened extra index for the Hz at max rsq (not FOI).
        % summaryRsqHzidMax_stride = nan(nSubs, 3);

        [rsqmax, hzmax_idstep] = max(rsq_step(1:hzlim1)); % due to nyquist.
        summaryRsqMax_step(ippant, ispeed) = rsqmax;
        summaryRsqHzidMax_step(ippant,ispeed) = Hzspace(hzmax_idstep);

        [rsqmax, hzmax_idstride] = max(rsq_stride(1:hzlim2));
        summaryRsqMax_stride(ippant, ispeed) = rsqmax;
        summaryRsqHzidMax_stride(ippant,ispeed) = Hzspace(hzmax_idstride);
        % add bound
        if ~isempty(nullCV_step)
            summaryNullCVMax_step(ippant, ispeed) = nullCV_step(3, hzmax_idstep);
        end
        if ~isempty(nullCV_stride)
            summaryNullCVMax_stride(ippant, ispeed) = nullCV_stride(3, hzmax_idstride);
        end

        %% Per-participant subplots -------------------------------------------
        if runParticipantPlots

            %% Column 1: step DV bar chart ------------------------------------
            ax1 = subplot(3, ncols, 1 + ncols*(ispeed-1));
            bar(ax1, xvec1, stepData, 'FaceColor', speedCols{ispeed}, 'FaceAlpha', 0.7);
            hold(ax1, 'on');



            ylabel(ax1, speedLabels{ispeed}, 'FontSize', fontsize, 'FontWeight', 'bold');
            title(ax1,  sprintf('Step %s (gc)', typeDV), 'FontSize', fontsize);
            xlabel(ax1, 'Gait cycle %', 'FontSize', fontsize);

            % add fit overlay?
            % Unbounded fourier1 best fit (a0 + a1·cos(w·x) + b1·sin(w·x))
            if overlayFIT
                [f, gof] = fit(xvec1', stepData', 'fourier1');
                h = plot(f, xvec1, stepData);
                h(2).LineWidth = 4;
                h(2).Color     = speedCols{ispeed};
                Hzapp      = xvec1(end) / (2*pi / f.w);
                legdetails = [sprintf('%.2f', Hzapp) ' Hz_{GC},  R^2 = ' ...
                    sprintf('%.2f', gof.rsquare)];
                legend(h(2), legdetails, 'FontSize', 15, ...
                    'AutoUpdate', 'off', 'Location', 'best');
                % repair labels.
                ylabel(ax1, speedLabels{ispeed}, 'FontSize', fontsize, 'FontWeight', 'bold');
                title(ax1,  sprintf('Step %s (gc)', typeDV), 'FontSize', fontsize);
                xlabel(ax1, 'Gait cycle %', 'FontSize', fontsize);

            end

            if ~isempty(dvYLim), ylim(ax1, dvYLim); end
            formatCountAxis(ax1, xvec1, fontsize);

            %% Column 2: step R²(Hz) ------------------------------------------
            ax2 = subplot(3, ncols, 2 + ncols*(ispeed-1));
            hold(ax2, 'on');
            plot(ax2, Hzspace, rsq_step, 'Color', speedCols{ispeed}, 'LineWidth', 2.5);

            if plotPermCV && ~isempty(nullCV_step)
                plot(ax2, Hzspace, nullCV_step(3, :), '--', ...
                    'Color', speedCols{ispeed}, 'LineWidth', 1.5);
            end


            xline(ax2, summaryRsqHzidMax_step(ippant,ispeed), ':', 'Color', [0.6 0.1 0.1], ...
                'LineWidth', 1.5, ...
                'Label', sprintf('%.2f cpstep', summaryRsqHzidMax_step(ippant,ispeed)), 'FontSize', 10);

            xlim(ax2, [0 nyqlim1]);
            ylim(ax2, [0 1]);
            xlabel(ax2, 'Frequency (cycles per step)', 'FontSize', fontsize);
            ylabel(ax2, 'R^2',                         'FontSize', fontsize);
            title(ax2, sprintf('Step R^2  (%.0f cps = %.2f)', ...
                summaryRsqHzidMax_step(ippant,ispeed), summaryRsqMax_step(ippant, ispeed)), ...
                'FontSize', fontsize);
            set(ax2, 'FontSize', fontsize);

            %% Column 3: stride DV bar chart ----------------------------------
            ax3 = subplot(3, ncols, 3 + ncols*(ispeed-1));
            bar(ax3, xvec2, strideData, 'FaceColor', speedCols{ispeed}, 'FaceAlpha', 0.7);
            hold(ax3, 'on');



            title(ax3,  sprintf('Stride %s (doubgc)', typeDV), 'FontSize', fontsize);
            xlabel(ax3, 'Gait cycle %', 'FontSize', fontsize);
            % Unbounded fourier1 best fit (a0 + a1·cos(w·x) + b1·sin(w·x))
            if overlayFIT
                [f, gof] = fit(xvec2', strideData', 'fourier1');
                h = plot(f, xvec2, strideData);
                h(2).LineWidth = 4;
                h(2).Color     = speedCols{ispeed};
                Hzapp      = xvec2(end) / (2*pi / f.w);
                legdetails = [sprintf('%.2f', Hzapp) ' Hz_{GC},  R^2 = ' ...
                    sprintf('%.2f', gof.rsquare)];
                legend(h(2), legdetails, 'FontSize', 15, ...
                    'AutoUpdate', 'off', 'Location', 'best');
                % repair labels.
                title(ax3,  sprintf('Stride %s (doubgc)', typeDV), 'FontSize', fontsize);
                xlabel(ax3, 'Gait cycle %', 'FontSize', fontsize);

            end
            if ~isempty(dvYLim), ylim(ax3, dvYLim); end
            formatCountAxis(ax3, xvec2, fontsize);

            %% Column 4: stride R²(Hz) ----------------------------------------
            ax4 = subplot(3, ncols, 4 + ncols*(ispeed-1));
            hold(ax4, 'on');
            plot(ax4, Hzspace, rsq_stride, 'Color', speedCols{ispeed}, 'LineWidth', 2.5);

            if plotPermCV && ~isempty(nullCV_stride)
                plot(ax4, Hzspace, nullCV_stride(3, :), '--', ...
                    'Color', speedCols{ispeed}, 'LineWidth', 1.5);
            end

            xline(ax4, summaryRsqHzidMax_stride(ippant,ispeed), ':', 'Color', [0.6 0.1 0.1], ...
                'LineWidth', 1.5, ...
                'Label', sprintf('%.2f cpstride', summaryRsqHzidMax_stride(ippant,ispeed)), 'FontSize', 10);

            xlim(ax4, [0 nyqlim2]);
            ylim(ax4, [0 1]);
            xlabel(ax4, 'Frequency (cycles per stride)', 'FontSize', fontsize);
            ylabel(ax4, 'R^2',                           'FontSize', fontsize);
            title(ax4, sprintf('Stride R^2  (%.0f cps = %.2f)', ...
                summaryRsqHzidMax_stride(ippant,ispeed), summaryRsqMax_stride(ippant, ispeed)), ...
                'FontSize', fontsize);
            set(ax4, 'FontSize', fontsize);

        end % runParticipantPlots

    end % ispeed

    %% Save per-participant figure ---------------------------------------------
    if runParticipantPlots
        sgtitle(['Participant ' num2str(ippant) ': ' subjIDs{ippant} ...
            '  |  ' typeDV ' — uniformity check  |  combined foot'], ...
            'FontSize', fontsize + 1, 'Interpreter', 'none');

        print(fullfile(inspectDir, ...
            ['Participant ' subjIDs{ippant} ' (' num2str(ippant) ') ' dvName ' uniformityCheck' appendp]), ...
            '-dpng');
        shg
        fprintf('Participant %d/%d (%s) saved.\n', ippant, nSubs, subjIDs{ippant});
    end

end % ippant


%% Summary figure: all participants × all speeds ==============================
%
%  Layout: 2 rows (step | stride) × 3 columns (Slow | Normal | Combined)
%  Participants are sorted by combined-speed stride R² (descending) so that
%  likely exclusions appear on the left — consistent with the original script.
%
if runSummaryPlot

    % Sort order: descending stride R² at combined speed
    [~, sortIdx] = sort(summaryRsq_stride(:, 3), 'descend');
    sortIdx = 1:size(summaryNullCV_stride, 1);   % or remove this line to use sorted order

    tickLabels = subjIDs(sortIdx);
    xpos       = 1:nSubs;

    figure(10); clf
    set(gcf, 'Color', 'w', 'Units', 'normalized', 'Position', [0.02 0.05 0.95 0.85]);

    for ispeed = 1:3

        rsq_step_sorted    = summaryRsq_step  (sortIdx, ispeed);
        rsq_stride_sorted  = summaryRsq_stride(sortIdx, ispeed);
        null_step_sorted   = summaryNullCV_step  (sortIdx, ispeed);
        null_stride_sorted = summaryNullCV_stride(sortIdx, ispeed);

        sig_step   = rsq_step_sorted   > null_step_sorted;
        sig_stride = rsq_stride_sorted > null_stride_sorted;

        %% Row 1: step R² at FOI — one column per speed -----------------------
        ax_top = subplot(2, 3, ispeed);
        bar(ax_top, xpos, rsq_step_sorted, 'FaceColor', speedCols{ispeed}, 'FaceAlpha', 0.8);
        hold(ax_top, 'on');

        if hasPFXNull
            for i = 1:nSubs
                if sig_step(i)
                    text(i, rsq_step_sorted(i) + 0.03, '*', ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

                    text(i, rsq_step_sorted(i) + 0.1, tickLabels{i}, ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

                end
            end
        end
        % add the PFX CV:

        plot(xpos, null_step_sorted, 'k:', 'LineWidth',2)


        set(ax_top, 'XTick', xpos, 'XTickLabel', tickLabels, ...
            'XTickLabelRotation', 90, 'FontSize', 11);
        ylabel(ax_top, sprintf('R^2 at %.0f cps', freqOfInterest_step), ...
            'FontSize', fontsize);
        title(ax_top, sprintf('Step — %s\n(%.0f cycle/step)', ...
            speedLabels{ispeed}, freqOfInterest_step), 'FontSize', fontsize);
        ylim(ax_top, [0 1]);
        box(ax_top, 'on');

        %% Row 2: stride R² at FOI — one column per speed --------------------
        ax_bot = subplot(2, 3, ispeed + 3);
        bar(ax_bot, xpos, rsq_stride_sorted, 'FaceColor', speedCols{ispeed}, 'FaceAlpha', 0.8);
        hold(ax_bot, 'on');

        if hasPFXNull
            for i = 1:nSubs
                if sig_stride(i)
                    text(i, rsq_stride_sorted(i) + 0.03, '*', ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k');

                    text(i, rsq_stride_sorted(i) + 0.1, tickLabels{i}, ...
                        'HorizontalAlignment', 'center', ...
                        'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
                end
            end
        end

        % add the PFX CV:
        plot(xpos, null_stride_sorted, 'k:', 'LineWidth',2)



        set(ax_bot, 'XTick', xpos, 'XTickLabel', tickLabels, ...
            'XTickLabelRotation', 90, 'FontSize', 11);
        ylabel(ax_bot, sprintf('R^2 at %.0f cps', freqOfInterest_stride), ...
            'FontSize', fontsize);
        title(ax_bot, sprintf('Stride — %s\n(%.0f cycles/stride)', ...
            speedLabels{ispeed}, freqOfInterest_stride), 'FontSize', fontsize);
        ylim(ax_bot, [0 1]);
        box(ax_bot, 'on');

    end % ispeed (summary columns)

    if hasPFXNull
        annotation('textbox', [0.72 0.01 0.25 0.03], 'String', ...
            '* p < .05 vs. null 95th-pct CV  (j9 full perm)', ...
            'EdgeColor', 'none', 'FontSize', 11, 'HorizontalAlignment', 'right');
    end

    sgtitle(sprintf( ...
        '%s — uniformity R^2  |  step: %.0f cps  |  stride: %.0f cps  |  sorted by combined stride R^2', ...
        typeDV, freqOfInterest_step, freqOfInterest_stride), ...
        'FontSize', fontsize + 1);

    print(fullfile(inspectDir, ...
        ['Summary_' dvName '_uniformityRanking_bySpeeds' appendp]), '-dpng');
    fprintf('\nSummary figure saved to: %s\n', inspectDir);

    %% now final summary plot, showing sig oscillations, mapped to frequencies.
    % so for this, use separate markers for sig/not [*, o],
    %y-axis can be Rsq, xAxis can be Hzspace.


    figure(11); clf
    set(gcf, 'Color', 'w', 'Units', 'normalized', 'Position', [0.02 0.05 0.95 0.85]);
    sortIdx = 1:size(summaryNullCV_stride, 1);   % or remove this line to use sorted order
    for ispeed = 1:3

        %now sort is needed, since we will place individually.
        rsq_step_sorted    = summaryRsqMax_step  (sortIdx, ispeed);
        rsq_stride_sorted  = summaryRsqMax_stride(sortIdx, ispeed);
        null_step_sorted   = summaryNullCVMax_step  (sortIdx, ispeed);
        null_stride_sorted = summaryNullCVMax_stride(sortIdx, ispeed);
        % also the CPS index:
        rsq_step_hz_sorted= summaryRsqHzidMax_step(sortIdx,ispeed);
        rsq_stride_hz_sorted= summaryRsqHzidMax_stride(sortIdx,ispeed);

        sig_step   = rsq_step_sorted   > null_step_sorted;
        sig_stride = rsq_stride_sorted > null_stride_sorted;

        %% Row 1: (LHS) max step R² at step CPS —  -----------------------
        ax_top = subplot(2, 6, ispeed);

        % sig
        plot(rsq_step_hz_sorted(sig_step), rsq_step_sorted(sig_step), 'r*');
        hold on;
        %non sig
        plot(rsq_step_hz_sorted(~sig_step), rsq_step_sorted(~sig_step), 'bo')

        % add text with subjID for debug:
        for ippant = 1:nSubs

            text(rsq_step_hz_sorted(ippant), rsq_step_sorted(ippant)+.02,subjIDs{ippant},...
                'fontsize', 7,...
                'HorizontalAlignment','center');
        end
        shg

        ylabel('Rsq')
        xlabel('Frequency (cycles per step)');
        title(ax_top, sprintf('Step — %s\n(Max Rsq in cps)', ...
            speedLabels{ispeed}), 'FontSize', fontsize);
        ylim(ax_top, [0 1]);
        box(ax_top, 'on');


%% also print list for easy ID: 
disp('*******')
disp(['Significant oscillations at step - ' speedLabels{ispeed}]);
disp('*******')
% copy paste from command window if needed.
%nan for non sig.
Sigoscillations= rsq_step_hz_sorted;
Sigoscillations(~sig_step)=nan;
tr= table((1:nSubs)',(subjIDs(:) ), rsq_step_hz_sorted, Sigoscillations);
disp(tr)

        %% add sorted list beneath for interpretation.
        subplot(2, 6, ispeed+6);
        % sort by max Rsq,
        [subRsq,sid] = sort(rsq_step_sorted, 'descend');
        bar(1:length(subRsq), subRsq);
        set(gca,'Xtick', 1:nSubs,  'XTickLabel', subjIDs(sid));
        % also add subjID, and HZ for interp.
        for ippant = 1:nSubs

            %were they sig?
            if sig_step(sid(ippant))
                scol='r';
            else
                scol='k';
            end
            % add approx cps
            text(ippant, subRsq(ippant)+.05, [sprintf('%.0f',rsq_step_hz_sorted(sid(ippant)))],...
                'HorizontalAlignment', 'center', 'color', scol)

            %add ID
            if mod(ippant,3)==1
                text(ippant, subRsq(ippant)+.1, [subjIDs(sid(ippant))],...
                    'HorizontalAlignment', 'center', 'Rotation', 90,'color', scol)
            end

        end
        ylim([0 1])
        %% Row 1: (RHS) max stride R² at stride CPS  --------------------
        ax_bot = subplot(2, 6, ispeed + 3);

        % sig
        plot(rsq_stride_hz_sorted(sig_stride), rsq_stride_sorted(sig_stride), 'r*');
        hold on;
        %non sig
        plot(rsq_stride_hz_sorted(~sig_stride), rsq_stride_sorted(~sig_stride), 'bo')

        % add text with subjID for debug:
        for ippant = 1:nSubs

            text(rsq_stride_hz_sorted(ippant), rsq_stride_sorted(ippant)+.02,subjIDs{ippant},...
                'fontsize', 7,...
                'HorizontalAlignment','center');
        end
        shg

        ylabel('Rsq')
        xlabel('Frequency (cycles per stride)');
        title(ax_bot, sprintf('Stride — %s\n(Max Rsq in cps)', ...
            speedLabels{ispeed}), 'FontSize', fontsize);
        ylim(ax_bot, [0 1]);
        box(ax_bot, 'on');

        %% also print list for easy ID: 
disp('*******')
disp(['Significant oscillations at STRIDE - ' speedLabels{ispeed}]);
disp('*******')
% copy paste from command window if needed.
%nan for non sig.
Sigoscillations= rsq_stride_hz_sorted;
Sigoscillations(~sig_stride)=nan;
tr= table((1:nSubs)',(subjIDs(:) ), rsq_stride_hz_sorted, Sigoscillations);
disp(tr)
        %% add sorted list beneath for interpretation.
        subplot(2, 6, ispeed+3 +6);
        % sort by max Rsq,
        [subRsq,sid] = sort(rsq_stride_sorted, 'descend');
        bar(1:length(subRsq), subRsq);
        set(gca,'Xtick', 1:nSubs,  'XTickLabel', subjIDs(sid));
        % also add subjID, and HZ for interp.
        for ippant = 1:nSubs

            %were they sig?
            if sig_stride(sid(ippant))
                scol='r';
            else
                scol='k';
            end
            % add approx cps
            text(ippant, subRsq(ippant)+.05, [sprintf('%.0f',rsq_stride_hz_sorted(sid(ippant)))],...
                'HorizontalAlignment', 'center', 'color', scol)

            %add ID
            if mod(ippant,3)==1
                text(ippant, subRsq(ippant)+.1, [subjIDs(sid(ippant))],...
                    'HorizontalAlignment', 'center', 'Rotation', 90,'color', scol)
            end

        end
        ylim([0 1])
    end % ispeed (summary columns)

    if hasPFXNull
        annotation('textbox', [0.72 0.01 0.25 0.03], 'String', ...
            '* p < .05 vs. null 95th-pct CV  (j9 full perm)', ...
            'EdgeColor', 'none', 'FontSize', 11, 'HorizontalAlignment', 'right');
    end

    sgtitle(sprintf( ...
        '%s — uniformity R^2  |  step: %.0f cps  |  stride: %.0f cps  ', ...
        typeDV, freqOfInterest_step, freqOfInterest_stride), ...
        'FontSize', fontsize + 1);

    print(fullfile(inspectDir, ...
        ['Summary_' dvName '_Significany by CPS and Speeds' appendp]), '-dpng');
    fprintf('\nSummary figure saved to: %s\n', inspectDir);

end % runSummaryPlot


% =========================================================================
%  LOCAL FUNCTION
% =========================================================================
function formatCountAxis(ax, xvec, fontsize)
midp = xvec(ceil(length(xvec) / 2));
set(ax, 'FontSize', fontsize, ...
    'XTick',      [xvec(1), midp, xvec(end)], ...
    'XTickLabel', {'0%', '50%', '100%'});
end

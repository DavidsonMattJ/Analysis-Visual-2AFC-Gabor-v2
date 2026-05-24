% j6_testfourier_Obs_NullGFX_vGaborv2.m
%
% Fourier-fit analysis: observed GFX data vs null distribution.
%
% For each DV × speed combination this script:
%   1. Fits a fourier1 model at each candidate frequency in Hzspace to the
%      observed group-mean binned time-series and stores R² at each frequency.
%   2. Optionally (testFull_perm = 1) fits the same model to every null
%      permutation and stores the [5th, 50th, 95th] percentile R² across
%      permutations as critical values for inference.
%
% Only stride-level data (doubgc, nGaits_toPlot = 2) is tested; two-cycle
% windows give better-constrained Fourier fits than single steps (gc).
%
% Inputs (from GFX/):
%   GFX_Data_inGaits[_new].mat          – observed binned SDT data  (j4b output)
%   GFX_Data_inGaits_null[_new].mat     – null shuffle distributions (j5 output)
%   GFX_Data_inGaits_FourierFits[_new].mat – prior results (appended if present)
%
% Saves to:
%   GFX/GFX_Data_inGaits_FourierFits[_new].mat
%     Variables: GFX_FourierNull, Hzspace
%
% Prerequisites: j4b_gaitPhaseResolved_vGaborv2, j5_createNull_testGFX_vGaborv2
% Runtime note: testFull_perm = 1 is slow; parfor runs over Hzspace (51 frequencies).

% -------------------------------------------------------------------------
%  Directories (guard: skip if homedir is already on the MATLAB path)
% -------------------------------------------------------------------------
if ~contains(path, homedir)
    set_myOnlineDirectories_vGaborv2;
end

% =========================================================================
%  CONTROL FLAG
% =========================================================================
% testFull_perm = 1 : also fits each null permutation → produces ShuffCV
%                     critical-value bands for inference against observed R².
% testFull_perm = 0 : observed-data frequency scan only (fast).
%                     ShuffCV output fields will be saved as all-NaN placeholders.
testFull_perm = 1; % will take approx 10 mins per DV

% -------------------------------------------------------------------------
%  File-name suffix (must match j4b / j5 output files)
% -------------------------------------------------------------------------
usenewTrialStride = 1;   % 1 = '_new' stride-epoch files; 0 = legacy files
if usenewTrialStride
    appendp = '_new';
else
    appendp = '';
end

%% Show participant list -------------------------------------------------------
cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs = length(pfols);
tr    = table((1:nsubs)', {pfols(:).name}');
disp(tr)

%% Load observed and null GFX data --------------------------------------------
cd([procdatadir filesep 'GFX'])
load('GFX_grandAvg_data.mat', 'GFX_headY'); % used for plotting
load(['GFX_Data_inGaits'      appendp]);   % observed binned SDT data (j4b output)
load(['GFX_Data_inGaits_null' appendp]);   % null shuffle data        (j5 output)

% Append to any previously saved Fourier-fit results, or start fresh.
try
    load(['GFX_Data_inGaits_FourierFits' appendp]);
catch
    GFX_FourierNull = [];   % no prior file – will be created on save
end

%% Number of permutations saved by j5 -----------------------------------------
% doubgc_Acc is [nPerm × 100]; first dimension gives the permutation count.
nPerm = size(GFX_TargPos_nullData(1,1,1).gc_Acc, 1);

%% Parallel pool ---------------------------------------------------------------
% Leave one core free for the OS; cap at 12 workers.
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Creating parallel pool...\n');
    parpool('local', min(12, feature('numcores') - 1));
else
    fprintf('Parallel pool already active: %d workers.\n', pool.NumWorkers);
end

% =========================================================================
%  DV SPECIFICATION
% =========================================================================
% Index:  1=Accuracy  2=RT  3=d'  4=criterion  5=HR  6=FA
%         7=counts (Target onset)  8=counts (Response onset)
%
% Note: index 8 is the only path that switches to the Response-onset dataset
% (GFX_RespPosData / GFX_RespPos_nullData). All other indices use Target onset.
% 'counts' appears at both 7 and 8 because the DV name is the same for both
% onset types; the onset dataset is selected by the testtype == 8 branch below.
%
% Extend the loop vector to add more DVs (e.g., add 7 for Target-onset counts).
% =========================================================================
testDVs = {'Accuracy','RT','dprime','crit','HR','FA','counts','counts'};

for testtype = [5];%,7]

    usebin = 1;   % 1 = binned data (recommended); 0 = unbinned (all 100 gait% points)

    % Default: align to Target onset ------------------------------------------
    dataIN    = GFX_TargPosData;
    datanull  = GFX_TargPos_nullData;
    typeOnset = 'Target';
    typeDV    = testDVs{testtype};

    % Exception: index 8 tests 'counts' aligned to Response onset ------------
    if testtype == 8
        dataIN    = GFX_RespPosData;
        datanull  = GFX_RespPos_nullData;
        typeOnset = 'Response';
        typeDV    = testDVs{testtype};
    end

    %% Build cfg ---------------------------------------------------------------
    cfg           = [];
    cfg.subjIDs   = subjIDs;
    cfg.type      = typeOnset;
    cfg.DV        = typeDV;   
    cfg.pidx1     = pidx1;
    cfg.pidx2     = pidx2;
    cfg.binDV     = 'Calc';         % 'Calc' = pooled-bin SDT; 'Av' = within-bin mean (Acc and RTs only).

    % Normalisation: off by default (RT data are already z-scored in j4b).
    % Set cfg.norm = 1 and choose normtype to re-enable.
    % cfg.ylims applies only when norm = 0.
    cfg.norm     = 1;
    cfg.normtype = 'relchange';   % options: 'absolute','relative','relchange','normchange','db'
    cfg.ylims    = [-.15  .15];

    % Participant exclusions applied before the group mean is computed.
    % Uncomment the minimal set (line 1) or use the extended list (line 2; currently active).
    
    cfg.omitPpants=[]; % none    
    % cfg.omitPpants= [27,36,32,8,31,3,19];% composite from poor gait
    % extraction
    %% Apply exclusions --------------------------------------------------------
    % Rows are removed from local copies only; source arrays are unaffected,
    % so the full dataset is available afresh on the next testtype iteration.
    if ~isempty(cfg.omitPpants)
        dataIN  (cfg.omitPpants, :, :) = [];
        datanull(cfg.omitPpants, :, :) = [];
    end

    % -------------------------------------------------------------------------
    %  Gait type: stride only (nGaits_toPlot = 2, i.e. doubgc)
    %  Two-cycle windows give better-constrained Fourier fits than single steps.
    %  To also test single-step fits, replace with:  for nGaits_toPlot = 1:2
    % -------------------------------------------------------------------------
    for nGaits_toPlot =1:2

    % Normalisation suffix appended to output field names ('norm' or empty).
    if cfg.norm == 1
        ntype = 'norm';
    else
        ntype = [];
    end

    %% Select struct field name for this DV ------------------------------------
    gaitfield = {'gc', 'doubgc'};
    binfield  = {'', '_binned'};

    if strcmp(cfg.DV, 'RT')
        usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts' cfg.binDV];
    elseif strcmp(cfg.DV, 'Accuracy')
        usefield = [gaitfield{nGaits_toPlot} '_binned_Acc' cfg.binDV];
    elseif strcmp(cfg.DV, 'crit')
        usefield = [gaitfield{nGaits_toPlot} '_binned_crit' cfg.binDV];
    else
        % Covers: dprime, HR, FA, counts (all share the same field-name pattern)
        usefield = [gaitfield{nGaits_toPlot} '_binned_' cfg.DV cfg.binDV];
    end

    % iLR = 3 → combined foot condition (no left/right split).
    % Consistent with j5, which only computes the stride-level combined condition.
    iLR = 3;

    %% Speed loop (1 = slow, 2 = normal, 3 = combined) ------------------------
    for iSpeed = 1:3

        %% Collate participant data: [nSubs × nBins] ---------------------------
        ppantData = [];
        shuffData = [];

        for isub = 1:size(dataIN, 1)
            ppantData(isub, :) = dataIN(isub, iSpeed, iLR).(usefield);
            if testFull_perm == 1
                % shuffData: [nSubs × nPerm × nBins]
                shuffData(isub, :, :) = datanull(isub, iSpeed, iLR).(usefield);
            end
        end % isub

        %% Optional normalisation (disabled when cfg.norm = 0) ----------------
        if cfg.norm == 1
            ppantData= applyNorm(ppantData, cfg.normtype);
        end

        %% Bin index and x-axis ------------------------------------------------
        if nGaits_toPlot == 1
            pidx    = cfg.pidx1;
            ftnames = {'LR', 'RL', 'combined'};   % retained for potential plotting use
        else
            pidx    = cfg.pidx2;
            ftnames = {'LRL', 'RLR', 'combined'};
        end

        % When using unbinned data, treat every sample point as its own bin.
        if usebin == 0
            pidx = 1:size(ppantData, 2);
        end

        % x-axis: approximate centre point of each bin
        mdiff = round(mean(diff(pidx)) / 2);
        xvec  = pidx(1:end-1) + mdiff;   % [1 × nBins]

        %% Group mean across participants (the series to be Fourier-fitted) ----
        gM = mean(ppantData, 1, 'omitnan');   % [1 × nBins]

        %% Fourier fit setup ---------------------------------------------------
        % Model: fourier1 = a0 + a1·cos(w·x) + b1·sin(w·x)
        % Strategy: pin w at each candidate frequency to obtain a frequency-
        % resolved R² profile, rather than allowing the optimiser to choose w freely.
        FitType    = 'fourier1';
        CoeffNames = coeffnames(fittype(FitType));

        % Bounds table: w bounds will be updated per frequency in the loop below.
        CoeffBounds = array2table( ...
            [-Inf(1, length(CoeffNames)); Inf(1, length(CoeffNames))], ...
            'RowNames',      ["lower bound", "upper bound"], ...
            'VariableNames', CoeffNames);

        % Candidate frequencies (Hz, relative to gait-cycle length): 51 values.
        Hzspace = 0.01:0.2:10;

        % Pre-allocate R² storage for the observed data sweep.
        fits_Rsquared_obsrvd = nan(1, length(Hzspace));

        % ShuffCV rows: [5th percentile; 50th (median); 95th percentile].
        % Remains all-NaN when testFull_perm = 0 (saved as placeholder for
        % downstream scripts that expect the field to exist).
        fits_Rsquared_shuffCV = nan(3, length(Hzspace));

        if testFull_perm == 1
            % Average shuffle data across subjects first: [nPerm × nBins]
            % shuffData is [nSubs × nPerm × nBins]; mean over dim 1.
            fits_Rsquared_shuff = nan(nPerm, length(Hzspace));
            meanShuff = squeeze(mean(shuffData, 1, 'omitnan'));   % [nPerm × nBins]
        end

        %% Pass 1: fit observed data at every candidate frequency --------------
        % FitOpts are built once here per frequency and stored so they can be
        % reused in the parallel null loop below, avoiding nPerm × nFreq
        % redundant fitoptions() calls.
        fitoptsperHz = cell(1, length(Hzspace));   % pre-allocated cell array

        for ifreq = 1:length(Hzspace)

            % Convert Hz (relative to gait cycle) → angular frequency w.
            % xvec(end) is treated as one full cycle period.
            testw = 2 * pi * Hzspace(ifreq) / xvec(end);

            % Pin w: set lower == upper bound to force the fit at this frequency.
            CoeffBounds.w(1) = testw;
            CoeffBounds.w(2) = testw;

            FitOpts = fitoptions('Method', 'NonlinearLeastSquares', ...
                'Lower', table2array(CoeffBounds(1, :)), ...
                'Upper', table2array(CoeffBounds(2, :)));

            [~, gof] = fit(xvec', gM', FitType, FitOpts);
            fits_Rsquared_obsrvd(1, ifreq) = gof.rsquare;

            fitoptsperHz{ifreq} = FitOpts;   % cache for null loop
        end % ifreq (observed)

        %% Pass 2 (optional): fit each null permutation ------------------------
        % parfor runs over Hzspace (51 frequencies) rather than nPerm (1000).
        % This avoids the overhead of re-launching the parallel pool 51 times;
        % each worker handles all permutations for one frequency serially.
        if testFull_perm == 1

            parfor ifreq = 1:length(Hzspace)
                usefit  = fitoptsperHz{ifreq};
                tmp_rsq = nan(nPerm, 1);

                for iperm = 1:nPerm
                    try
                        [~, gof] = fit(xvec', squeeze(meanShuff(iperm, :))', ...
                            FitType, usefit);
                        tmp_rsq(iperm) = gof.rsquare;
                    catch
                        % Most likely a NaN in the null mean series; leave as NaN.
                    end
                end % iperm

                % Sliced assignment: each ifreq writes to a separate column (parfor-safe).
                fits_Rsquared_shuff(:, ifreq) = tmp_rsq;
            end % parfor ifreq

            % Compute percentile critical values across permutations at each frequency.
            for ifreq = 1:length(Hzspace)
                fits_Rsquared_shuffCV(:, ifreq) = ...
                    quantile(fits_Rsquared_shuff(:, ifreq), [.05, .5, .95]);
            end

            fprintf('  Null fits complete  —  %s | %s | speed %d\n', ...
                typeOnset, typeDV, iSpeed);
        end % testFull_perm

        %% Store results in output struct ---------------------------------------
        obsField = [cfg.type 'Ons_' usefield '_fitsRsq_Obs'     ntype];
        cvField  = [cfg.type 'Ons_' usefield '_fitsRsq_ShuffCV' ntype];
        GFX_FourierNull(iSpeed).(obsField) = fits_Rsquared_obsrvd;
        GFX_FourierNull(iSpeed).(cvField)  = fits_Rsquared_shuffCV;

        fprintf('  Done  —  gait %d | %s | %s | speed %d/3\n', ...
            nGaits_toPlot, typeOnset, typeDV, iSpeed);

    end % iSpeed
    end % nGaits
end % testtype

%% Save (append Hzspace alongside fit results) --------------------------------
fprintf('\nSaving Fourier fit results...\n');
cd([procdatadir filesep 'GFX']);
try save(['GFX_Data_inGaits_FourierFits' appendp], ...
    'GFX_FourierNull', 'Hzspace', '-append');
catch
    save(['GFX_Data_inGaits_FourierFits' appendp], ...
    'GFX_FourierNull', 'Hzspace');
end

fprintf('Saved: GFX_Data_inGaits_FourierFits%s.mat\n', appendp);

%%
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
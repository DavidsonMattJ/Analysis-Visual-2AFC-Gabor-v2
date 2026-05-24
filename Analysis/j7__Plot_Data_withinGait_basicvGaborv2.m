% j7__Plot_Data_withinGait_basicvGaborv2.m
%
% Master plotting script: gait-phase-resolved SDT time-series (within-gait).
% Dispatches to one or more job sections depending on the toggle flags below.
%
% Job sections
% ------------
%   plotforInspection     – per-participant QC plots for exclusion decisions
%                           (convenience shortcut; canonical call is in run_pipeline.m)
%   plotBasic_superscript – main multipurpose gait-phase plot (group level or PFX, all
%                           foot × speed combinations via plot_GaitresultsBinned_2speed_vAUD)
%   plotforACNS_2024      – poster/conference layout: 4 DVs in a row at one speed,
%                           plus an auditory vs visual sensitivity overlay panel
%   plot_GrandAverageResults – grand-average bar/scatter plots across walk speeds
%                              (convenience shortcut; canonical call is in run_pipeline.m)
%
% Inputs (loaded from GFX/):
%   GFX_Data_inGaits[_new].mat          – binned SDT data         (j4b output)
%   GFX_Data_inGaits_FourierFits[_new].mat – Fourier fit results  (j6/j9 output)
%
% Prerequisites: j4b_gaitPhaseResolved_vAUD
% Optional:      j6_testfourier_Obs_NullGFX_vAud, j9__testfourier_Obs_NullPFX_vAUD

% -------------------------------------------------------------------------
%  Directories (guard: skip if homedir is already on the MATLAB path)
% -------------------------------------------------------------------------
% if ~contains(path, homedir)
    set_myOnlineDirectories_vGaborv2;
% end

% =========================================================================
%  TOGGLE FLAGS
% =========================================================================
% Set a field to 1 to run that job section, 0 to skip.
%
%   plotforInspection     – calls plotPFX_forInspection; also run directly
%                           in Stage 5 of run_pipeline.m (set 0 to avoid duplication)
%   plotBasic_superscript – main within-gait grouped plot
%   plotforACNS_2024      – conference poster layout (requires setmydirs_detectvGabor)
%   plot_GrandAverageResults – calls job_MSplot2_grandaverage; also run directly
%                              in Stage 5 of run_pipeline.m (set 0 to avoid duplication)
% =========================================================================
jobs = [];
jobs.plotforInspection      = 0;
jobs.plotBasic_superscript  = 1;
jobs.plotforACNS_2024       = 0;
jobs.plot_GrandAverageResults = 0;

% -------------------------------------------------------------------------
%  Options
% -------------------------------------------------------------------------
% omitbadPpants: apply the standard exclusion list before any plotting.
% Exclusions are applied to local working copies so source arrays stay intact
% and each job section receives identically filtered data.
omitbadPpants = 1;

% File-name suffix — guard so the pipeline value is not overwritten if set upstream.
if ~exist('usenewTrialStride', 'var')
    usenewTrialStride = 1;
end
if usenewTrialStride
    appendp = '_new';
else
    appendp = '';
end

%% Show participant list -------------------------------------------------------
cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
tr    = table((1:length(pfols))', {pfols(:).name}');
disp(tr)

%% Load GFX data --------------------------------------------------------------
fprintf('Loading GFX data...\n');
cd([procdatadir filesep 'GFX'])
load(['GFX_Data_inGaits' appendp]);   % binned SDT data (j4b output)
load('GFX_grandAvg_data', 'GFX_headY') ; % for plot overlays (j4a output).

% Fourier fit results saved by j6 (GFX) and j9 (PFX) into the same file.
% Optional — downstream jobs degrade gracefully if absent.
try
    load(['GFX_Data_inGaits_FourierFits' appendp]);
catch
    warning(['GFX_Data_inGaits_FourierFits' appendp '.mat not found. ' ...
        'Fourier fit overlays (plotFITforced, plotShuff) will be unavailable. ' ...
        'Run j6 (GFX) and/or j9 (PFX) to generate this file.']);
end

%% Apply participant exclusions -----------------------------------------------
% Exclusions are applied to local working copies of each array.
% GFX_FourierNull is a group-level result — exclusion was applied at compute time
% (j6), so it is not modified here.
% PFX_FourierNull is participant-level (j9 output) — rows are removed to match.

%criteria summarised in participantInclusions.xls

omitPpants=[];

if omitbadPpants
    GFX_TargPosData_plot = GFX_TargPosData;   % working copies
    GFX_RespPosData_plot = GFX_RespPosData;
    GFX_headY_plot       = GFX_headY;
    subjIDs_plot         = subjIDs;

    GFX_TargPosData_plot(omitPpants, :, :) = [];
    GFX_RespPosData_plot(omitPpants, :, :) = [];
    GFX_headY_plot      (omitPpants, :, :) = [];
    subjIDs_plot        (omitPpants)        = [];

    % PFX_FourierNull is only available after j9 has run.
    if exist('PFX_FourierNull', 'var')
        PFX_FourierNull_plot = PFX_FourierNull;
        PFX_FourierNull_plot(omitPpants, :) = [];
    else
        PFX_FourierNull_plot = [];
    end

    fprintf('Exclusions applied: %d participants removed, %d remaining.\n', ...
        numel(omitPpants), size(GFX_TargPosData_plot, 1));
else
    % No exclusions — working copies point to the original arrays.
    GFX_TargPosData_plot = GFX_TargPosData;
    GFX_RespPosData_plot = GFX_RespPosData;
    GFX_headY_plot       = GFX_headY;
    subjIDs_plot         = subjIDs;
    PFX_FourierNull_plot = [];
    if exist('PFX_FourierNull', 'var')
        PFX_FourierNull_plot = PFX_FourierNull;
    end
end


% =========================================================================
%% JOB: plotBasic_superscript
% =========================================================================
% Main within-gait grouped plot via plot_GaitresultsBinned_2speed.
% Output: 3 × 3 figure grid — columns = speeds, rows = foot conditions.
% Covers all foot (L, R, combined) × speed (slow, normal, combined)
% combinations for the chosen DV.
%
% cfg.DV options: 'Accuracy' | 'RT' | 'dprime' | 'crit' | 'counts' | 'FA' | 'HR'
% cfg.binDV:      'Calc' (pooled-bin SDT recalculation) | 'Av' (within-bin mean)
%
% Note: PFX_FourierNull is passed as [] — participant-level Fourier inference
% is handled by j10 (plot_Data_PFX_FourierFits). GFX_FourierNull provides
% the group-level fit overlays here.
if jobs.plotBasic_superscript

    % clf
    cfg = [];
    cfg.plotlevel    = 'GFX';          % 'GFX' = group; 'PFX' = per participant
    cfg.usebin       = 1;              % 1 = binned; 0 = all 100 gait% points
    cfg.type         = 'Target';       % 'Target' | 'Response' (Response not fully supported)
    cfg.DV           = 'HR';           % see options above
    cfg.binDV        = 'Calc';         % 'Calc' = pooled-bin SDT; 'Av' = within-bin mean
    cfg.plotCOL      = 'b';
    cfg.fitCOL       = 'b';
    cfg.yyaxis       = 'left';
    cfg.ispeed       = 1:3;            % speeds to include: 1=slow, 2=normal, 3=combined
    cfg.nGaitstoplot = 1:2;              % 1 = single step (gc), 2 = stride (doubgc)
    cfg.datadir      = datadir;
    cfg.HeadData     = GFX_headY_plot;
    cfg.pidx1        = pidx1;
    cfg.pidx2        = pidx2;
    cfg.subjIDs      = subjIDs_plot;
    cfg.figdir       = figdir;
    cfg.normON       = 1;
    cfg.normtype     = 'relchange';    % 'relative' | 'relchange' | 'normchange' | 'absolute'
    cfg.plotHead     = 0;              % 1 = overlay head-position trace
    cfg.plotFIT      = 1;             % 1 = overlay best (unbounded) Fourier fit
    cfg.plotFITforced = 1;            % 1 = overlay frequency-resolved R² fits (j6/j9)
    cfg.plotShuff    = 0;             % 1 = overlay null CV bands (requires j6 full perm)
    cfg.omitPpants   = [];            % exclusions already applied above; leave empty here
    cfg.plotGaitGrid = 0;
    cfg.plotnewLayout = 0;            % 0 = all foot × gait combinations; 1 = focused subset

    % PFX_FourierNull passed as [] — participant-level Fourier overlays are
    % intentionally omitted here; use j10 for PFX Fourier figures.
    plot_GaitresultsBinned_2speed_vGaborv2(GFX_TargPosData_plot, GFX_FourierNull, cfg);
%%
    plotdebug_ppantoverlay;
end % plotBasic_superscript


% =========================================================================
%% JOB: plotforACNS_2024
% =========================================================================
% Conference poster layout (ACNS 2024 / EPC 2025).
% Part A: 4 DVs in a single row at one fixed speed (ispeed = 2, Normal).
% Part B: Auditory vs visual sensitivity overlay — reloads visual (Gabor)
%         data alongside auditory data, then calls the overlaid plot function.
%
% NOTE: Part B calls setmydirs_detectvGabor, which overwrites procdatadir,
% datadir, and related path variables. Auditory directories are restored
% immediately afterwards via set_myOnlineDirectories_AUD. All visual data
% are stored in local variables (GFX_TargPos_Data_gabor etc.) before the
% directory switch so no auditory workspace variables are corrupted.
if jobs.plotforACNS_2024

    %% Part A: 4 DVs in a row (auditory, one speed) ---------------------------
    clf
    cfg = [];
    cfg.plotlevel   = 'GFX';
    cfg.usebin      = 1;
    cfg.type        = 'Target';                            % 'Target' | 'Response'
    cfg.DV_all      = {'Accuracy', 'HR', 'dprime', 'crit'};
    cfg.plotCOL_all = {'b', 'r', 'm', [.5 .5 .9]};
    cfg.fitCOL_all  = {'b', 'r', 'm', [.5 .5 .9]};
    cfg.yyaxis_all  = {'left', 'right', 'left', 'right'};
    cfg.ispeed      = 2;                                   % Normal walking speed only
    cfg.datadir     = datadir;
    cfg.HeadData    = GFX_headY_plot;
    cfg.pidx1       = pidx1;
    cfg.pidx2       = pidx2;
    cfg.subjIDs     = subjIDs_plot;
    cfg.figdir      = figdir;
    cfg.plotShuff   = 0;
    cfg.add_cpsFit  = 0;
    cfg.plotHead    = 0;
    cfg.plotFIT     = 1;

    for iDV = 1:4
        cfg.DV      = cfg.DV_all{iDV};
        cfg.plotCOL = cfg.plotCOL_all{iDV};
        cfg.fitCOL  = cfg.fitCOL_all{iDV};
        cfg.yyaxis  = cfg.yyaxis_all{iDV};
        cfg.iDV     = iDV;

        if strcmp(cfg.type, 'Target')
            % Three args: data, GFX null fits, cfg
            plot_GaitresultsBinned_2speed_ACNS(GFX_TargPosData_plot, GFX_FourierNull, cfg);
        else
            % Response onset: GFX_FourierNull was computed on Target onset;
            % pass [] to avoid mismatched field access.
            plot_GaitresultsBinned_2speed_ACNS(GFX_RespPosData_plot, [], cfg);
        end
    end % iDV

    %% Part B: Auditory vs visual d' and criterion overlay --------------------
    % Step 1: load visual (Gabor) data into local variables.
    % setmydirs_detectvGabor overwrites procdatadir/datadir — auditory dirs
    % are restored immediately after the load.
    setmydirs_detectvGabor;
    cd([procdatadir filesep 'GFX']);
    load('GFX_Data_inGaits.mat', 'GFX_TargPosData', 'pidx2');
    GFX_TargPos_Data_gabor = GFX_TargPosData;
    GFXNull_gabor          = GFX_FourierNull;
    pidx2_gab              = pidx2;

    % Step 2: restore auditory directories and load auditory data.
    set_myOnlineDirectories_AUD;
    cd([procdatadir filesep 'GFX']);
    load('GFX_Data_inGaits.mat', 'GFX_TargPosData', 'pidx1', 'pidx2', 'subjIDs');
    GFX_TargPos_Data_aud = GFX_TargPosData;
    GFXNull_aud          = GFX_FourierNull;
    pidx2_aud            = pidx2;

    % Step 3: configure and plot overlay (visual dashed, auditory solid).
    clf
    cfg = [];
    cfg.plotlevel     = 'GFX';
    cfg.usebin        = 1;
    cfg.type          = 'Target';
    cfg.DV_all        = {'dprime', 'dprime', 'crit', 'crit'};
    cfg.DV_modalityall = {'visual', 'auditory', 'visual', 'auditory'};
    cfg.plotCOL_all   = {'m', 'm', [.5 .5 .9], [.5 .5 .9]};
    cfg.fitCOL_all    = {'m', 'm', [.5 .5 .9], [.5 .5 .9]};
    cfg.yyaxis_all    = {'left', 'left', 'left', 'left'};
    cfg.FITstyle      = {'-', ':', '-', ':'};   % solid = visual, dashed = auditory
    cfg.ispeed        = 3;                       % combined speed
    cfg.datadir       = datadir;
    cfg.HeadData      = GFX_headY_plot;          % auditory head data (post-exclusion)
    cfg.pidx1         = pidx1;
    cfg.pidx2         = pidx2;
    cfg.subjIDs       = subjIDs_plot;
    cfg.figdir        = figdir;
    cfg.plotShuff     = 0;
    cfg.add_cpsFit    = 0;
    cfg.plotHead      = 0;
    cfg.plotFIT       = 1;

    for iDV = 1:4
        cfg.DV         = cfg.DV_all{iDV};
        cfg.plotCOL    = cfg.plotCOL_all{iDV};
        cfg.fitCOL     = cfg.fitCOL_all{iDV};
        cfg.yyaxis     = cfg.yyaxis_all{iDV};
        cfg.DV_modality = cfg.DV_modalityall{iDV};
        cfg.iDV        = iDV;

        if iDV == 1 || iDV == 3   % visual panels
            GFX_T        = GFX_TargPos_Data_gabor;
            GFX_N        = GFXNull_gabor;
            cfg.pidx2    = pidx2_gab;
        else                       % auditory panels (iDV == 2 || iDV == 4)
            GFX_T        = GFX_TargPos_Data_aud;
            GFX_N        = GFXNull_aud;
            cfg.pidx2    = pidx2_aud;
        end

        plot_GaitresultsBinned_2speed_ACNS_overlaid(GFX_T, GFX_N, cfg);
    end % iDV

end % plotforACNS_2024


% =========================================================================
%% JOB: plot_GrandAverageResults
% =========================================================================
% Grand-average bar/scatter plots comparing walk speeds (not split by stride).
% This is a convenience shortcut — the canonical call is made directly in
% Stage 5 of run_pipeline.m. Set to 0 when running the full pipeline to
% avoid generating duplicates.
if jobs.plot_GrandAverageResults
    job_MSplot2_grandaverage;
end

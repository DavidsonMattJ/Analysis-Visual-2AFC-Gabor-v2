% job_MSplot3_intraStride
%
% -------------------------------------------------------------------------
%  Options
% -------------------------------------------------------------------------
% omitbadPpants: apply the standard exclusion list before any plotting.
% Exclusions are applied to local working copies so source arrays stay intact
% and each job section receives identically filtered data.




% set_myOnlineDirectories_AUD; omitbadPpants = 1;
set_myOnlineDirectories_vGaborv2; omitBadPpants=0; 

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

% omitPpants=[22,6,15,5,28,2, 4,16,     17,20];
omitPpants= [2 6 20 21 22 27 28 30 32 38];

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
%% JOB: MS figure
% =========================================================================
% Main within-gait grouped plot via plot_GaitresultsBinned_2speed_vAUD_overlay.
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
% if jobs.plotBasic_superscript

    % clf
    cfg = [];
    cfg.plotlevel    = 'GFX';          % 'GFX' = group; 'PFX' = per participant
    cfg.usebin       = 1;              % 1 = binned; 0 = all 100 gait% points
    cfg.type         = 'Target';       % 'Target' | 'Response' (Response not fully supported)
    cfg.DV           = 'dprime';           % see options above
    cfg.binDV        = 'Calc';         % 'Calc' = pooled-bin SDT; 'Av' = within-bin mean
    cfg.plotCOL      = 'k';            %for linking bins.
    cfg.fitCOL       = 'k';            %override with speedCols.
    cfg.yyaxis       = 'left';
    cfg.ispeed       = 1:2;            % speeds to include: 1=slow, 2=normal, 3=combined
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
    cfg.plotFITforced = 0;            % 1 = overlay frequency-resolved R² fits (j6/j9)
    cfg.plotShuff    = 0;             % 1 = overlay null CV bands (requires j6 full perm)
    cfg.omitPpants   = [];            % exclusions already applied above; leave empty here
    
    % PFX_FourierNull passed as [] — participant-level Fourier overlays are
    % intentionally omitted here; use j10 for PFX Fourier figures.
    
    plot_GaitresultsBinned_2speed_vGabor_overlay(GFX_TargPosData_plot, GFX_FourierNull, cfg);

% end % plotBasic_superscript


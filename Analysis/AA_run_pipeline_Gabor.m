%% run_pipeline.m
% Master pipeline script for the Gabor Discrimination (while walking) analysis.
%
% PIPELINE STAGES
% ---------------
%   STAGE 0  -- Setup & data repair (run once before everything else)
%   STAGE 1  -- Import & gait detection   (j0, j1)
%   STAGE 2  -- Gait cycle annotation     (j2, j2B, j2C, j3)
%   STAGE 3  -- Group-level crunching     (j4a, j4b)
%   STAGE 4  -- Statistical testing       (j5, j6, j9)
%   STAGE 5  -- Figures                   (j7, j8, j10, grand-average plots)
%
% Toggle the jobs struct fields below to run selected stages.
% Set a field to 1 to run, 0 to skip.
%
% =========================================================================

clear jobs

% ---- Stage 0: directories (always runs) ---------------------------------
set_myOnlineDirectories_vGaborv2;

% =========================================================================
%  TOGGLE STAGES HERE
% =========================================================================
jobs.stage1_import        = 0;   % j0, j1
jobs.stage2_gaitAnnot     = 0;   % j2, j2B, j2C, j3
jobs.stage3_crunch        = 0;   % j4a, j4b
jobs.stage4_stats         = 0;   % j5, j6, j9
jobs.stage5_figures       = 0;   % j7, j8, j10 + grand-average plots
% =========================================================================


%% STAGE 1 – Import raw VR data and detect gait peaks
% j0: reads raw CSVs → HeadPos + trial_summaryTable (per participant)
% j1: detects gait cycle peaks/troughs in head position
if jobs.stage1_import
    fprintf('\n=== STAGE 1: Import & gait detection ===\n');
    
    j0_ImportVRData;
    % from here, we can inspect individual data for staircase fails, etc.
    
    plotPFX_forInspection_vGaborv2;
    %
    j1_findpeaksinHeadPos;
end


%% STAGE 2 – Annotate gait cycles and epoch time-series
% j2:  appends gait-cycle percentage and SDT labels (Hit/Miss/FA/CR)
% j2B: classifies left/right foot placement per gait cycle
% j2C: extends to double gait cycles (strides: LRL / RLR)
% j3:  epochs and resamples head-position time-series onto gait grid
if jobs.stage2_gaitAnnot
    fprintf('\n=== STAGE 2: Gait cycle annotation ===\n');
    j2_appendgaitcycleData_wSDT_vGaborv2; 
    j2B_appendgaitcycledata_wLR_vGaborv2;
    j2C_append_doubGCdata_gaborv2;
    j3_epochGait_timeseries_SlowNormal_vGaborv2;
end


%% STAGE 3 – Concatenate participants and bin into group-level structure
% j4a: grand-average head timeseries + SDT metrics per walk speed
%      → saves GFX/GFX_grandAvg_data.mat
% j4b: gait-phase-resolved SDT, binned by step/stride × speed × foot
%      → saves GFX/GFX_Data_inGaits_new.mat
if jobs.stage3_crunch
    fprintf('\n=== STAGE 3: Group crunching ===\n');
    j4a_grandAverage_vGaborv2;
    j4b_gaitPhaseResolved_vGaborv2;
end


%% STAGE 4 – Statistical testing (Fourier + permutation null)
% j5:  builds null distributions via gait-position shuffling (1000 perms)
% j6:  Fourier analysis – group level (GFX), observed vs null
% j9:  Fourier analysis – participant level (PFX), observed vs null
%      (slow; uses parfor – ensure parallel toolbox is available)
if jobs.stage4_stats
    fprintf('\n=== STAGE 4: Statistical testing ===\n');
    j5_createNull_testGFX_vGaborv2;
    j6_testfourier_Obs_NullGFX_vGaborv2; % full perm waiting.

    % j9__testfourier_Obs_NullPFX_vAUD;   % slow – comment out if not needed
end


%% STAGE 5 – Figures (inspection, exclusion)
% j7:  within-gait plots (binned DV timecourses with Fourier fits)
% j8:  Fourier-fit summary plots – group level
% j10: Fourier-fit + phase plots – participant level
% grand-average bar/box plots (two style options)

if jobs.stage5_figures
    fprintf('\n=== STAGE 5: Figures ===\n');

    % Inspection plot – useful before finalising figure scripts:
    plotPFX_forInspection; % staircases and Acc, HR, FA.
    % plotPFX_exclusionCheck_targetCounts% gaitresolved
    % plotPFX_exclusionCheck_byDV; % gaitresolved
    % 
    % % Main figure scripts:
    % j7__Plot_Data_withinGait_basicvAUD;
    % j8__Plot_Data_FourierFits_GFX_vAUD;
    % j10__Plot_Data_PFX_FourierFits_phase_prev_vAud;
    % 
    % % Grand-average speed comparison figures (choose one or both):
    % job_MSplot2_grandaverage;   % box-and-scatter style
    % plot_MSfigure2;             % raincloud style
end


%% STAGE 6 - Figures (MS ready).
if jobs.stage6_MSFigures
    fprintf('\n=== STAGE 5: Figures ===\n');
    
    % job_MSplot2_grandAverageSummary; % shows step dur, Acc, RT (box and scatter).
    % 
    % job_MSplot3_intraStride_v1; % shows step and stride for HR, cps for HR and targ counts overlaid

end


fprintf('\n=== Pipeline complete ===\n');

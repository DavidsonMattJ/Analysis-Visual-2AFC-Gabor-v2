% j3_epochGait_timeseries_SlowNormal_vGaborv2

% using the trough information in the summary table (or HeadPos), we will
% now epoch the raw, resampled, and normalized versions of the gait
% time-series (for later plots).

% Auditory version - Two walking speeds

%add paths if new session/ running in isolation.
if ~contains(path, homedir)
    set_myOnlineDirectories_vGaborv2;
end

cd(procdatadir)
%% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs = length(pfols);
tr = table((1:length(pfols))', {pfols(:).name}');
disp(tr)
%%

resampSize = 100;

for ippant = 1:nsubs
    cd(procdatadir)
    load(pfols(ippant).name, ...
        'HeadPos', 'trial_summaryTable', 'subjID');
    savename = pfols(ippant).name;

    disp(['Preparing j3 gait time series data... ' savename]);

    %% find all gait onsets (in samples) and get bad trial list once
    gOnsets   = find(~isnan(trial_summaryTable.trgO_gStart_samp));
    ngaits    = length(gOnsets);
    badtrials = rejTrials_detectvGaborv2(subjID);

    %% Pre-allocate output arrays with NaN (no zero-to-NaN cleanup needed)
    gait_ts_raw        = nan(ngaits, 500);
    gait_ts_resamp     = nan(ngaits, resampSize);
    doubgait_ts_raw    = nan(ngaits, 500);
    doubgait_ts_resamp = nan(ngaits, 200);

    %% Pre-allocate metadata table
    gait_ts_gData = table( ...
        nan(ngaits,1), nan(ngaits,1), nan(ngaits,1), nan(ngaits,1), ...
        repmat({'L'}, ngaits, 1), nan(ngaits,1), nan(ngaits,1), nan(ngaits,1), ...
        'VariableNames', {'trialallocation','walkSpeed','SDTcat','gaitIdx', ...
                          'gaitFeet','gaitStart','gaitSamps','gaitDuration'});

    %% Single pass over all gait onsets
    for igait = 1:ngaits

        rowIndx = gOnsets(igait);
        itrial  = trial_summaryTable.trial(rowIndx);

        if ismember(itrial, badtrials)
            continue
        end

        %--- Detrend once, reuse for both single and double gait ---
        headY = detrend(HeadPos(itrial).Y);

        %--- Single gait cycle ---
        gaitStFin = [trial_summaryTable.trgO_gStart_samp(rowIndx), ...
                     trial_summaryTable.trgO_gFin_samp(rowIndx)];
        gaitDur   = trial_summaryTable.trgO_gFin_sec(rowIndx) - ...
                    trial_summaryTable.trgO_gStart_sec(rowIndx);
        gFt       = trial_summaryTable.trgO_gFoot(rowIndx);

        rawHead    = headY(gaitStFin(1):gaitStFin(2));
        resampHead = imresize(rawHead', [1, resampSize]);
        gIdx       = find(HeadPos(itrial).Y_gait_troughs == gaitStFin(1));

        gait_ts_raw(igait, 1:length(rawHead)) = rawHead;
        gait_ts_resamp(igait, :)              = resampHead;

        gait_ts_gData.trialallocation(igait) = itrial;
        gait_ts_gData.SDTcat(igait)          = trial_summaryTable.SDTcat(rowIndx);
        gait_ts_gData.walkSpeed(igait)       = trial_summaryTable.walkSpeed(rowIndx);
        gait_ts_gData.gaitDuration(igait)    = gaitDur;
        gait_ts_gData.gaitFeet(igait)        = gFt;
        gait_ts_gData.gaitIdx(igait)         = gIdx;
        gait_ts_gData.gaitStart(igait)       = gaitStFin(1);
        gait_ts_gData.gaitSamps(igait)       = length(rawHead);

        %--- Double gait cycle (steps n and n+1) ---
        allGs = HeadPos(itrial).Y_gait_troughs;
        try
            doubStFin = allGs(gIdx) : allGs(gIdx+2);
        catch
            doubStFin = allGs(gIdx-1) : allGs(gIdx+1);
        end

        rawDoub    = headY(doubStFin);
        resampDoub = imresize(rawDoub', [1, 200]);

        doubgait_ts_raw(igait, 1:length(rawDoub)) = rawDoub;
        doubgait_ts_resamp(igait, :)              = resampDoub;

    end % igait

    disp(['saving gait and double gait time series data for ' subjID])
    cd(procdatadir);
    save(savename, 'gait_ts_raw', 'gait_ts_resamp', 'gait_ts_gData', ...
        'doubgait_ts_raw', 'doubgait_ts_resamp', '-append');

end % per ppant

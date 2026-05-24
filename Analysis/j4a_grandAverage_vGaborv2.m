% j4a_grandAverage_vGaborv2
%
% Computes per-participant grand average metrics:
%   - Mean head position timeseries per walk speed (single and double GC)
%   - Average gait duration
%   - SDT metrics (accuracy, RT, HR, FAR, d-prime, criterion) per walk speed
%
% Output saved to: GFX/GFX_grandAvg_data.mat
%   GFX_headY(ippant, iSpeed).gc / .doubgc  -- mean resampled timeseries
%     iSpeed: 1=slow, 2=normal, 3=combined
%   GFX_grandAvg(ippant, iSpeed).grand_*    -- scalar SDT metrics
%   avWalkParams(ippant)                    -- mean stride duration (L+R step)
%   subjIDs                                 -- cell array of participant IDs
%
% Run before j4b_gaitPhaseResolved_vAUD.

if ~contains(path, homedir)
    set_myOnlineDirectories_vGaborv2
end

cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs = length(pfols);
tr = table((1:length(pfols))', {pfols(:).name}');
disp(tr)

%% Initialise outputs
GFX_headY    = [];
GFX_grandAvg = [];
subjIDs      = cell(1, nsubs);
avStepDuration = nan(nsubs,2); % two speeds.

for ippant = 1:nsubs
    cd(procdatadir)
    load(pfols(ippant).name, 'subjID', 'trial_summaryTable', ...
        'gait_ts_gData', 'gait_ts_resamp', 'doubgait_ts_resamp');

    subjIDs{ippant} = subjID;
    disp(['Grand average: ' subjID]);

    %% Head position timeseries averages
    slowsteps   = gait_ts_gData.walkSpeed == 1;
    normalsteps = gait_ts_gData.walkSpeed == 2;

    GFX_headY(ippant, 1).gc     = mean(gait_ts_resamp(slowsteps,   :), 1, 'omitnan');
    GFX_headY(ippant, 2).gc     = mean(gait_ts_resamp(normalsteps, :), 1, 'omitnan');
    GFX_headY(ippant, 3).gc     = mean(gait_ts_resamp,              1, 'omitnan');

    GFX_headY(ippant, 1).doubgc = mean(doubgait_ts_resamp(slowsteps,   :), 1, 'omitnan');
    GFX_headY(ippant, 2).doubgc = mean(doubgait_ts_resamp(normalsteps, :), 1, 'omitnan');
    GFX_headY(ippant, 3).doubgc = mean(doubgait_ts_resamp,              1, 'omitnan');

    % Average stride duration: sum of mean L-step and mean R-step durations
    Ltrials_ts = strcmp(gait_ts_gData.gaitFeet, 'LR');
    Rtrials_ts = strcmp(gait_ts_gData.gaitFeet, 'RL');
    
    avStepDuration(ippant,1) = mean(gait_ts_gData.gaitDuration(slowsteps), 'omitnan'); % note this is only on trials with targets.
    avStepDuration(ippant,2) = mean(gait_ts_gData.gaitDuration(normalsteps), 'omitnan');
                           

    %% Grand average SDT per walk speed
    slowspeedTrials   = find(trial_summaryTable.walkSpeed == 1);
    normalspeedTrials = find(trial_summaryTable.walkSpeed == 2);
    allspeedTrials    = sort(unique([slowspeedTrials; normalspeedTrials]));
    speedstoIndex     = {slowspeedTrials, normalspeedTrials, allspeedTrials};

    sigPres = find(trial_summaryTable.signalPresent == 1);
    sigAbs  = find(trial_summaryTable.signalPresent == 0);

    for iSpeed = 1:3
        usespd = speedstoIndex{iSpeed};

        grand_acc = mean(trial_summaryTable.targCor(usespd),  'omitnan');
        grand_rt  = mean(trial_summaryTable.reactionTime(usespd),  'omitnan');

        % Hit rate
        trialselect = intersect(usespd, sigPres);
        allcats     = trial_summaryTable.SDTcat(trialselect);
        HITn        = sum(allcats == 1);
        grand_HR    = HITn / length(trialselect);

        % FA rate
        trialselect = intersect(usespd, sigAbs);
        allcats     = trial_summaryTable.SDTcat(trialselect);
        FAn         = sum(allcats == 3);
        grand_FAR   = FAn / length(trialselect);

        % Macmillan & Kaplan (1985) correction: replace 0 with 1/(2N), 1 with 1-1/(2N)
        Ntarg = sum(~isnan(trial_summaryTable.targOnset(usespd)));
        Nresp = sum(~isnan(trial_summaryTable.reactionTime(usespd)));

        if grand_HR  == 1, grand_HR  = 1 - 1/(2*Ntarg); end
        if grand_HR  == 0, grand_HR  =     1/(2*Ntarg); end
        if grand_FAR == 1, grand_FAR = 1 - 1/(2*Nresp); end
        if grand_FAR == 0, grand_FAR =     1/(2*Nresp); end

        grand_dprime = norminv(grand_HR) - norminv(grand_FAR);
        grand_crit   = -0.5 * (norminv(grand_HR) + norminv(grand_FAR));

        if isinf(grand_dprime) || isinf(grand_crit)
            warning('j4a: infinite d-prime or criterion for %s speed %d', subjID, iSpeed);
        end

        GFX_grandAvg(ippant, iSpeed).grand_Acc    = grand_acc;
        GFX_grandAvg(ippant, iSpeed).grand_rt     = grand_rt;
        GFX_grandAvg(ippant, iSpeed).grand_HR     = grand_HR;
        GFX_grandAvg(ippant, iSpeed).grand_FAR    = grand_FAR;
        GFX_grandAvg(ippant, iSpeed).grand_dprime = grand_dprime;
        GFX_grandAvg(ippant, iSpeed).grand_crit   = grand_crit;

    end % iSpeed

end % ippant

%% Save
cd([procdatadir filesep 'GFX']);
disp('Saving grand average data (GFX_grandAvg_data.mat)');
save('GFX_grandAvg_data', 'GFX_headY', 'GFX_grandAvg', 'avStepDuration', 'subjIDs');

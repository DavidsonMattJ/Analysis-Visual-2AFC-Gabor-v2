% j5_createNull_testGFX_vGaborv2
%
% Generates null distributions for the gait-phase-resolved SDT metrics by
% repeatedly shuffling the gait-position assignments across trials.
%
% Mirrors j4b_gaitPhaseResolved_vGaborv2 exactly for data selection, binning,
% and per-bin SDT recalculation. The only difference is that onset labels
% are randomly permuted (trial-level) before computing DVs at each gait%.
%
% Counts additionally use a bin-level shuffle of the empirical histogram.
%
% nPerm = 1000 shuffles per participant × condition combination.
% Expected runtime: ~120 s per participant (per DV , 3 speeds, step+stride).
% Expected runtime: ~40 s per participant (stride only - skips L/R assignment).
%
% Saves to: GFX/GFX_Data_inGaits_null_new.mat
%
% Run j4b_gaitPhaseResolved_vGaborv2 before this script.

if ~contains(path, homedir)
    set_myOnlineDirectories_vGaborv2;
end

cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs = length(pfols);
tr = table((1:length(pfols))', {pfols(:).name}');
disp(tr)

%% Load binning parameters from observed data (must match j4b exactly)
cd([procdatadir filesep 'GFX'])
load('GFX_Data_inGaits_new', 'pidx1', 'pidx2', 'gaittypes');

% Append to any previously saved null results, or start fresh.
try
    load(['GFX_Data_inGaits_null_new']);
catch
    GFX_TargPos_nullData = [];
    GFX_RespPos_nullData = []; % no prior file – will be created on save
end

cd(procdatadir)

nPerm = 1000;
nbins = length(pidx1) - 1;  % = 10

%%
for ippant = 1:nsubs
    tic
    cd(procdatadir)
    load(pfols(ippant).name, 'trial_summaryTable', 'subjID');
    ppantData = trial_summaryTable;

    disp(['Constructing null distributions: ' subjID]);

    %% Speed and foot indices — identical to j4b
    slowspeedTrials   = find(ppantData.walkSpeed == 1);
    normalspeedTrials = find(ppantData.walkSpeed == 2);
    allspeedTrials    = sort(unique([slowspeedTrials; normalspeedTrials]));
    speedstoIndex     = {slowspeedTrials, normalspeedTrials, allspeedTrials};

    Ltrials     = find(strcmp(ppantData.trgO_gFoot, 'LR'));
    Rtrials     = find(strcmp(ppantData.trgO_gFoot, 'RL'));
    LRAlltrials = sort(unique([Ltrials; Rtrials]));
    LRtoIndex   = {Ltrials, Rtrials, LRAlltrials};

    for nGait = 1:2   % 1 = single step,  2 = stride

        pidx = {pidx1, pidx2};
        pidx = pidx{nGait};

        % Stride: combined foot only (no L/R split)
        if nGait == 1
            
            LRloop = 1:3;
        else
            LRloop = 3;
        end

        for iSpeed = 1:3
            for iLR = LRloop

                %% Select rows — identical to j4b
                usespd = speedstoIndex{iSpeed};

                if nGait == 1
                    uset = intersect(LRtoIndex{iLR}, usespd);
                    tmp_targOnsets = ppantData.trgO_gPcnt(uset);
                    tmp_respOnsets = ppantData.respO_gPcnt(uset);
                else
                    strideRows     = find(~isnan(ppantData.trgO_gPcnt_strideinTrial));
                    uset           = intersect(strideRows, usespd);
                    tmp_targOnsets = ppantData.trgO_gPcnt_strideinTrial(uset);
                    tmp_respOnsets = ppantData.respO_gPcnt_strideinTrial(uset);
                end

                tmp_Correct = ppantData.targCor(uset);
                tmp_clickRT = ppantData.reactionTime(uset);
                tmp_SDTcat  = ppantData.SDTcat(uset);
                nTrials     = numel(uset);

                %% Pre-allocate null arrays: [nPerm × 100]
                [targOns_Acc_shuff,  respOns_Acc_shuff, ...
                 targOns_RTs_shuff,  respOns_RTs_shuff] = deal(nan(nPerm, 100));

                % H / M / FA / CR counts per gait% per permutation
                [targOns_SDTcat_shuff, respOns_SDTcat_shuff] = deal(nan(nPerm, 100, 4));

                %% Permutation loop: shuffle onset assignment, recompute DVs
                for iperm = 1:nPerm

                    % Trial-level randomisation: shuffle which trial's onset
                    % lands at which gait%, keeping DVs paired with trials.
                    permTargOnsets = tmp_targOnsets(randperm(nTrials));
                    permRespOnsets = tmp_respOnsets(randperm(nTrials));

                    for ip = 1:100
                        tIdx = find(permTargOnsets == ip);
                        rIdx = find(permRespOnsets == ip);

                        targOns_Acc_shuff(iperm, ip) = mean(tmp_Correct(tIdx), 'omitnan');
                        respOns_Acc_shuff(iperm, ip) = mean(tmp_Correct(rIdx), 'omitnan');
                        targOns_RTs_shuff(iperm, ip) = mean(tmp_clickRT(tIdx), 'omitnan');
                        respOns_RTs_shuff(iperm, ip) = mean(tmp_clickRT(rIdx), 'omitnan');

                        targOns_SDTcat_shuff(iperm, ip, 1) = sum(tmp_SDTcat(tIdx) == 1); % H
                        targOns_SDTcat_shuff(iperm, ip, 2) = sum(tmp_SDTcat(tIdx) == 2); % M
                        targOns_SDTcat_shuff(iperm, ip, 3) = sum(tmp_SDTcat(tIdx) == 3); % FA
                        targOns_SDTcat_shuff(iperm, ip, 4) = sum(tmp_SDTcat(tIdx) == 4); % CR

                        respOns_SDTcat_shuff(iperm, ip, 1) = sum(tmp_SDTcat(rIdx) == 1);
                        respOns_SDTcat_shuff(iperm, ip, 2) = sum(tmp_SDTcat(rIdx) == 2);
                        respOns_SDTcat_shuff(iperm, ip, 3) = sum(tmp_SDTcat(rIdx) == 3);
                        respOns_SDTcat_shuff(iperm, ip, 4) = sum(tmp_SDTcat(rIdx) == 4);
                    end

                end % iperm

                %% Count null: bin-level shuffle of the empirical histogram
                % Complements the trial-level null above: randomises which
                % bin gets which empirical count, preserving the marginal
                % distribution of event density across the gait cycle.
                edges = 0.5 : 1 : (pidx(end) + 0.5);   % bins centred on 1..pidx(end)
                [targOns_Counts_shuff, respOns_Counts_shuff] = deal(nan(nPerm, pidx(end)));

                for itype = 1:2
                    if itype == 1
                        useposdata = tmp_targOnsets;
                    else
                        useposdata = tmp_respOnsets;
                    end

                    real_counts = histcounts(useposdata, edges);   % 1 × pidx(end)
                    nbins_count = length(real_counts);
                    outgoingnull = nan(nPerm, nbins_count);

                    for iperm = 1:nPerm
                        outgoingnull(iperm, :) = real_counts(randperm(nbins_count));
                    end

                    if itype == 1
                        targOns_Counts_shuff = outgoingnull;
                    else
                        respOns_Counts_shuff = outgoingnull;
                    end
                end

                %% Bin shuffled metrics: [nPerm × nbins]
                % Boundary rule: pidx(i):pidx(i+1)-1  (matches j4b exactly)
                [targOns_Acc_binAv_shuff,    respOns_Acc_binAv_shuff, ...
                 targOns_RT_binAv_shuff,      respOns_RT_binAv_shuff,  ...
                 targOns_Counts_binAv_shuff, respOns_Counts_binAv_shuff,...
                  targOns_Counts_binCalc_shuff, respOns_Counts_binCalc_shuff,... 
                 targOns_RT_binCalc_shuff,respOns_RT_binCalc_shuff]  = deal(nan(nPerm, nbins));

                for ibin = 1:nbins
                    idx = pidx(ibin):(pidx(ibin+1)-1);
                    targOns_Acc_binAv_shuff(:, ibin)     = mean(targOns_Acc_shuff(:, idx),    2, 'omitnan');
                    respOns_Acc_binAv_shuff(:, ibin)     = mean(respOns_Acc_shuff(:, idx),    2, 'omitnan');
                    targOns_RT_binAv_shuff(:, ibin)      = mean(targOns_RTs_shuff(:, idx),    2, 'omitnan');
                    respOns_RT_binAv_shuff(:, ibin)      = mean(respOns_RTs_shuff(:, idx),    2, 'omitnan');
                    targOns_Counts_binAv_shuff(:, ibin)  = mean(targOns_Counts_shuff(:, idx),  2,'omitnan');
                    respOns_Counts_binAv_shuff(:, ibin)  = mean(respOns_Counts_shuff(:, idx),  2, 'omitnan');


                    targOns_RT_binCalc_shuff(:,ibin) = mean(targOns_RTs_shuff(:,idx),2, 'omitnan'); % matches the above for later scripts.
                    targOns_Counts_binCalc_shuff(:,ibin) = sum(targOns_Counts_shuff(:,idx),2, 'omitnan');
                    respOns_RT_binCalc_shuff(:,ibin) = mean(respOns_RTs_shuff(:,idx), 2,   'omitnan'); % matches the above for later scripts.
                    respOns_Counts_binCalc_shuff(:,ibin) = sum(respOns_Counts_shuff(:,idx),2,'omitnan');


                end

                %% Per-bin SDT recalculation from pooled shuffled counts
                % Matches j4b: Macmillan & Kaplan (1985) correction using
                % mean bin count as N (consistent across all bins).
                Ntargstridebin = nan(1, nbins);
                for ibin = 1:nbins
                    idx = pidx(ibin):(pidx(ibin+1)-1);
                    binCounts = sum(targOns_SDTcat_shuff(:, idx, :), [2 3]); % [nPerm×1×1]
                    Ntargstridebin(ibin) = mean(binCounts(:), 'omitnan');
                end
                Ncorr = mean(Ntargstridebin, 'omitnan');

                onsetData = {targOns_SDTcat_shuff, respOns_SDTcat_shuff};

                for itype = 1:2

                    usedata = onsetData{itype};
                    [dprime_bin_shuff, criterion_bin_shuff, ...
                     HitRate_bin_shuff, FARate_bin_shuff, Acc_bin_shuff] = deal(nan(nPerm, nbins));

                    for ibin = 1:nbins
                        idx = pidx(ibin):(pidx(ibin+1)-1);

                        HITn_bin  = sum(usedata(:, idx, 1), 2);  % [nPerm × 1]
                        MISSn_bin = sum(usedata(:, idx, 2), 2);
                        FAn_bin   = sum(usedata(:, idx, 3), 2);
                        CRn_bin   = sum(usedata(:, idx, 4), 2);

                        Hratebin  = HITn_bin  ./ (HITn_bin  + MISSn_bin);
                        FAratebin = FAn_bin   ./ (FAn_bin   + CRn_bin);

                        % Macmillan & Kaplan correction (vectorised over nPerm)
                        Hratebin(Hratebin  == 1) = 1 - 1/(2*Ncorr);
                        Hratebin(Hratebin  == 0) =     1/(2*Ncorr);
                        FAratebin(FAratebin == 1) = 1 - 1/(2*Ncorr);
                        FAratebin(FAratebin == 0) =     1/(2*Ncorr);

                        HitRate_bin_shuff(:, ibin)   = Hratebin;
                        FARate_bin_shuff(:, ibin)    = FAratebin;
                        dprime_bin_shuff(:, ibin)    = norminv(Hratebin)  - norminv(FAratebin);
                        criterion_bin_shuff(:, ibin) = -0.5*(norminv(Hratebin) + norminv(FAratebin));
                        Acc_bin_shuff(:, ibin)       = (HITn_bin + CRn_bin) ./ ...
                                                       (HITn_bin + CRn_bin + FAn_bin + MISSn_bin);
                    end % ibin

                    if itype == 1
                        targOns_Acc_binCalc_shuff       = Acc_bin_shuff;
                        targOns_dprime_binCalc_shuff    = dprime_bin_shuff;
                        targOns_criterion_binCalc_shuff = criterion_bin_shuff;
                        targOns_HitRate_binCalc_shuff   = HitRate_bin_shuff;
                        targOns_FARate_binCalc_shuff    = FARate_bin_shuff;
                    else
                        respOns_Acc_binCalc_shuff       = Acc_bin_shuff;
                        respOns_dprime_binCalc_shuff    = dprime_bin_shuff;
                        respOns_criterion_binCalc_shuff = criterion_bin_shuff;
                        respOns_HitRate_binCalc_shuff   = HitRate_bin_shuff;
                        respOns_FARate_binCalc_shuff    = FARate_bin_shuff;
                    end

                end % itype

                %% Store in output structs (field names mirror j4b)
                gn = {'gc_', 'doubgc_'};
                gn = gn{nGait};

                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'counts'])            = targOns_Counts_shuff;
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'Acc'])               = targOns_Acc_shuff;
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'rts'])               = targOns_RTs_shuff;
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'binned_countsAv'])   = targOns_Counts_binAv_shuff;
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'binned_countsCalc']) = targOns_Counts_binCalc_shuff;               
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'binned_AccAv'])      = targOns_Acc_binAv_shuff;
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'binned_AccCalc'])      = targOns_Acc_binCalc_shuff;
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'binned_rtsAv'])      = targOns_RT_binAv_shuff;
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'binned_rtsCalc'])      = targOns_RT_binCalc_shuff;
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'binned_dprimeCalc']) = targOns_dprime_binCalc_shuff;
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'binned_critCalc'])   = targOns_criterion_binCalc_shuff;
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'binned_HRCalc'])     = targOns_HitRate_binCalc_shuff;
                GFX_TargPos_nullData(ippant, iSpeed, iLR).([gn 'binned_FACalc'])     = targOns_FARate_binCalc_shuff;

                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'counts'])            = respOns_Counts_shuff;
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'Acc'])               = respOns_Acc_shuff;
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'rts'])               = respOns_RTs_shuff;
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'binned_countsAv']) = respOns_Counts_binAv_shuff;
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'binned_countsCalc']) = respOns_Counts_binCalc_shuff;
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'binned_AccAv'])      = respOns_Acc_binAv_shuff;
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'binned_AccCalc'])      = respOns_Acc_binCalc_shuff;
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'binned_rtsAv'])      = respOns_RT_binAv_shuff;
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'binned_rtsCalc'])      = respOns_RT_binCalc_shuff;                
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'binned_dprimeCalc']) = respOns_dprime_binCalc_shuff;
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'binned_critCalc'])   = respOns_criterion_binCalc_shuff;
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'binned_HRCalc'])     = respOns_HitRate_binCalc_shuff;
                GFX_RespPos_nullData(ippant, iSpeed, iLR).([gn 'binned_FACalc'])     = respOns_FARate_binCalc_shuff;

            end % iLR
            fprintf('  Speed %d/3 complete\n', iSpeed);
        end % iSpeed
    end % nGait

    fprintf('Participant %d (%s) done in %.1f s\n', ippant, subjID, toc);

end % ippant

%% Save
cd([procdatadir filesep 'GFX']);
disp('Saving null distribution data (GFX_Data_inGaits_null_new.mat)');
save('GFX_Data_inGaits_null_new', ...
    'GFX_TargPos_nullData', 'GFX_RespPos_nullData', ...
    'nPerm', 'pidx1', 'pidx2', 'gaittypes');

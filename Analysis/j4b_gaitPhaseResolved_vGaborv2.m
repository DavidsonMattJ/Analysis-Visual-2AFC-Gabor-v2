% j4b_gaitPhaseResolved_vGaborv2
%
% Computes gait-phase-resolved SDT metrics per participant, binned across
% 10 equal intervals of the step (1-100%) and stride (1-100%) cycle.
%
% Output struct dimensions:
%   GFX_TargPosData(ippant, iSpeed, iLR)
%   GFX_RespPosData(ippant, iSpeed, iLR)
%     iSpeed: 1=slow, 2=normal, 3=combined
%     iLR:    1=LR foot (single step only)
%             2=RL foot (single step only)
%             3=combined feet / all strides
%
%   Stride (doubgc) fields are only populated at iLR=3 (no L/R distinction
%   for strides after j2C simplification).
%
% Binning boundary rule: bin i covers pidx(i) : pidx(i+1)-1.
% Gait% index 100 is excluded from all bins (effectively never observed).
%
% Saved to: GFX/GFX_Data_inGaits_new.mat
%
% Run j4a_grandAverage_vGaborv2 first.

if ~contains(path, homedir)
    set_myOnlineDirectories_vGaborv2;
end

cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs = length(pfols);
tr = table((1:length(pfols))', {pfols(:).name}');
disp(tr)

%% Binning setup
% 10 bins spanning gait% 1-100.
% Using fewer bins than the 1-100 range due to trial counts per walk speed.
pidx1  = ceil(linspace(1, 101, 11));  % single step: 10 bins
pidx2  = ceil(linspace(1, 101, 11));  % stride:      10 bins
nbins  = length(pidx1) - 1;           % = 10

gaittypes = {'single gait', 'double gait'};

GFX_TargPosData = [];
GFX_RespPosData = [];
subjIDs         = cell(1, nsubs);

for ippant = 1:nsubs
    cd(procdatadir)
    load(pfols(ippant).name, 'subjID', 'trial_summaryTable');

    subjIDs{ippant} = [subjID '_' num2str(ippant)];
    disp(['Gait-phase analysis: ' subjID]);

    ppantData = trial_summaryTable;

    %% Speed and foot indices
    slowspeedTrials   = find(ppantData.walkSpeed == 1);
    normalspeedTrials = find(ppantData.walkSpeed == 2);
    allspeedTrials    = sort(unique([slowspeedTrials; normalspeedTrials]));
    speedstoIndex     = {slowspeedTrials, normalspeedTrials, allspeedTrials};

    % L/R foot classification (single step only; from j2B)
    Ltrials     = find(strcmp(ppantData.trgO_gFoot, 'LR'));
    Rtrials     = find(strcmp(ppantData.trgO_gFoot, 'RL'));
    LRAlltrials = sort(unique([Ltrials; Rtrials]));
    LRtoIndex   = {Ltrials, Rtrials, LRAlltrials};

    for nGait = 1:2   % 1 = single step,  2 = stride

        pidx = {pidx1, pidx2};
        pidx = pidx{nGait};

        % Stride has no L/R split — only compute combined (iLR=3)
        if nGait == 1
            LRloop = 1:3;
        else
            LRloop = 3;
        end

        for iSpeed = 1:3
            for iLR = LRloop

                %% Select rows for this speed × foot × gait-type combination
                usespd = speedstoIndex{iSpeed};

                if nGait == 1
                    uset = intersect(LRtoIndex{iLR}, usespd);
                    tmp_targOnsets = ppantData.trgO_gPcnt(uset);
                    tmp_respOnsets = ppantData.respO_gPcnt(uset);
                else
                    strideRows = find(~isnan(ppantData.trgO_gPcnt_strideinTrial));
                    uset       = intersect(strideRows, usespd);
                    tmp_targOnsets = ppantData.trgO_gPcnt_strideinTrial(uset);
                    tmp_respOnsets = ppantData.respO_gPcnt_strideinTrial(uset);
                end

                tmp_Correct = ppantData.targCor(uset);
                tmp_clickRT = ppantData.reactionTime(uset);
                tmp_SDTcat  = ppantData.SDTcat(uset);

                %% Per gait% (1-100): compute acc, RT, and SDT counts
                [targOns_Acc, respOns_Acc, ...
                 targOns_RTs, respOns_RTs, ...
                 targOns_Counts, respOns_Counts] = deal(nan(1, 100));

                % H / M / FA / CR counts at each gait% index
                [targOns_SDTcat, respOns_SDTcat] = deal(nan(1, 100, 4));

                for ip = 1:100
                    tIdx = find(tmp_targOnsets == ip);
                    rIdx = find(tmp_respOnsets == ip);

                    targOns_Acc(ip)    = mean(tmp_Correct(tIdx), 'omitnan');
                    respOns_Acc(ip)    = mean(tmp_Correct(rIdx), 'omitnan');
                    targOns_RTs(ip)    = mean(tmp_clickRT(tIdx), 'omitnan');
                    respOns_RTs(ip)    = mean(tmp_clickRT(rIdx), 'omitnan');
                    targOns_Counts(ip) = numel(tIdx);
                    respOns_Counts(ip) = numel(rIdx);

                    targOns_SDTcat(1, ip, 1) = sum(tmp_SDTcat(tIdx) == 1); % Hit
                    targOns_SDTcat(1, ip, 2) = sum(tmp_SDTcat(tIdx) == 2); % Miss
                    targOns_SDTcat(1, ip, 3) = sum(tmp_SDTcat(tIdx) == 3); % FA
                    targOns_SDTcat(1, ip, 4) = sum(tmp_SDTcat(tIdx) == 4); % CR

                    respOns_SDTcat(1, ip, 1) = sum(tmp_SDTcat(rIdx) == 1);
                    respOns_SDTcat(1, ip, 2) = sum(tmp_SDTcat(rIdx) == 2);
                    respOns_SDTcat(1, ip, 3) = sum(tmp_SDTcat(rIdx) == 3);
                    respOns_SDTcat(1, ip, 4) = sum(tmp_SDTcat(rIdx) == 4);
                end

                %% Bin simple metrics: mean acc/RT, summed counts
                % Boundary rule: bin i = pidx(i) : pidx(i+1)-1 (no double-counting).
                [targOns_Acc_binAv,    respOns_Acc_binAv, ...
                 targOns_RT_binAv,     respOns_RT_binAv,  ...
                 targOns_Counts_binAv, respOns_Counts_binAv,...
                 targOns_Counts_binCalc, respOns_Counts_binCalc,... 
                 targOns_RT_binCalc,respOns_RT_binCalc] = deal(nan(1, nbins));

                for ibin = 1:nbins
                    
                    idx = pidx(ibin):(pidx(ibin+1)-1);
                    

                    targOns_Acc_binAv(ibin)     = mean(targOns_Acc(idx),    'omitnan');
                    respOns_Acc_binAv(ibin)     = mean(respOns_Acc(idx),    'omitnan');
                    targOns_RT_binAv(ibin)      = mean(targOns_RTs(idx),    'omitnan');
                    respOns_RT_binAv(ibin)      = mean(respOns_RTs(idx),    'omitnan');
                    targOns_Counts_binAv(ibin) = mean(targOns_Counts(idx),  'omitnan');
                    respOns_Counts_binAv(ibin) = mean(respOns_Counts(idx),  'omitnan');

                    % can also compute (non-average) calculations at the bin level,
                    % for non SDT metrics:
                    targOns_RT_binCalc(ibin) = mean(targOns_RTs(idx),    'omitnan'); % matches the above for later scripts.
                    targOns_Counts_binCalc(ibin) = sum(targOns_Counts(idx), 'omitnan');
                    respOns_RT_binCalc(ibin) = mean(respOns_RTs(idx),    'omitnan'); % matches the above for later scripts.
                    respOns_Counts_binCalc(ibin) = sum(respOns_Counts(idx),'omitnan');

                    %accuracy at whole bin level calculated below based on
                    %all trial types (SDT cats) within.
                end

                %% Per-bin SDT recalculation (pooled raw counts within bin)
                % d-prime and criterion are computed from pooled H/M/FA/CR
                % counts in each bin, not by averaging per-index values.
                % Macmillan & Kaplan (1985): replace HR/FAR of 0 with 1/(2N)
                % and 1 with 1-1/(2N), where N = mean bin count across stride.

                % Mean bin count used as N for the M&K correction
                Ntargstridebin = nan(1, nbins);
                for ibin = 1:nbins
                    
                    idx = pidx(ibin):(pidx(ibin+1)-1);
                    
                    Ntargstridebin(ibin) = sum(targOns_SDTcat(1, idx, :), 'all');
                end
                Ncorr = mean(Ntargstridebin, 'omitnan');

                onsetData = {targOns_SDTcat, respOns_SDTcat};

                for itype = 1:2   % 1 = target onset,  2 = response onset

                    usedata = onsetData{itype};
                    [dprime_bin, criterion_bin, HitRate_bin, FARate_bin, Acc_bin] = ...
                        deal(nan(1, nbins));

                    for ibin = 1:nbins
                        idx = pidx(ibin):(pidx(ibin+1)-1);

                        HITn_bin  = sum(usedata(1, idx, 1), 2);
                        MISSn_bin = sum(usedata(1, idx, 2), 2);
                        FAn_bin   = sum(usedata(1, idx, 3), 2);
                        CRn_bin   = sum(usedata(1, idx, 4), 2);

                        Hratebin  = HITn_bin  / (HITn_bin  + MISSn_bin);
                        FAratebin = FAn_bin   / (FAn_bin   + CRn_bin);

                        % Macmillan & Kaplan correction
                        if Hratebin  == 1, Hratebin  = 1 - 1/(2*Ncorr); end
                        if Hratebin  == 0, Hratebin  =     1/(2*Ncorr); end
                        if FAratebin == 1, FAratebin = 1 - 1/(2*Ncorr); end
                        if FAratebin == 0, FAratebin =     1/(2*Ncorr); end

                        HitRate_bin(ibin)   = Hratebin;
                        FARate_bin(ibin)    = FAratebin;
                        dprime_bin(ibin)    = norminv(Hratebin)  - norminv(FAratebin);
                        criterion_bin(ibin) = -0.5*(norminv(Hratebin) + norminv(FAratebin));
                        Acc_bin(ibin)       = (HITn_bin + CRn_bin) / ...
                                              (HITn_bin + CRn_bin + FAn_bin + MISSn_bin);
                        % can also st
                        if isinf(dprime_bin(ibin))
                            warning('j4b: inf d-prime: %s spd%d LR%d nGait%d bin%d', ...
                                subjID, iSpeed, iLR, nGait, ibin);
                        end

                    end % ibin

                    if itype == 1
                        targOns_Accuracy_binCalc  = Acc_bin;
                        targOns_dprime_binCalc    = dprime_bin;
                        targOns_criterion_binCalc = criterion_bin;
                        targOns_HitRate_binCalc   = HitRate_bin;
                        targOns_FARate_binCalc    = FARate_bin;
                    else
                        respOns_Accuracy_binCalc  = Acc_bin;
                        respOns_dprime_binCalc    = dprime_bin;
                        respOns_criterion_binCalc = criterion_bin;
                        respOns_HitRate_binCalc   = HitRate_bin;
                        respOns_FARate_binCalc    = FARate_bin;
                    end

                end % itype

                %% Store in output structs
                gn = {'gc_', 'doubgc_'};
                gn = gn{nGait};

                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'counts'])            = targOns_Counts;
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'Acc'])               = targOns_Acc;
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'rts'])               = targOns_RTs;
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'binned_countsAv'])   = targOns_Counts_binAv;
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'binned_countsCalc']) = targOns_Counts_binCalc;
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'binned_AccAv'])      = targOns_Acc_binAv;
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'binned_AccCalc'])    = targOns_Accuracy_binCalc;
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'binned_rtsAv'])      = targOns_RT_binAv;
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'binned_rtsCalc'])      = targOns_RT_binCalc;                
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'binned_dprimeCalc']) = targOns_dprime_binCalc;
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'binned_critCalc'])   = targOns_criterion_binCalc;
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'binned_HRCalc'])     = targOns_HitRate_binCalc;
                GFX_TargPosData(ippant, iSpeed, iLR).([gn 'binned_FACalc'])     = targOns_FARate_binCalc;

                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'counts'])            = respOns_Counts;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'Acc'])               = respOns_Acc;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'rts'])               = respOns_RTs;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'binned_countsAv']) = respOns_Counts_binAv;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'binned_countsCalc']) = respOns_Counts_binCalc;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'binned_AccAv'])      = respOns_Acc_binAv;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'binned_rtsAv'])      = respOns_RT_binAv;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'binned_rtsCalc'])      = respOns_RT_binCalc;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'binned_AccCalc'])    = respOns_Accuracy_binCalc;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'binned_dprimeCalc']) = respOns_dprime_binCalc;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'binned_critCalc'])   = respOns_criterion_binCalc;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'binned_HRCalc'])     = respOns_HitRate_binCalc;
                GFX_RespPosData(ippant, iSpeed, iLR).([gn 'binned_FACalc'])     = respOns_FARate_binCalc;

            end % iLR
        end % iSpeed
    end % nGait

end % ippant

%% Save
% Filename kept as GFX_Data_inGaits_new for downstream compatibility.
% Note: grand average metrics are in GFX_grandAvg_data.mat (from j4a).
cd([procdatadir filesep 'GFX']);
disp('Saving gait-phase resolved data (GFX_Data_inGaits_new.mat)');
save('GFX_Data_inGaits_new', ...
    'GFX_TargPosData', 'GFX_RespPosData', ...
    'subjIDs', 'pidx1', 'pidx2', 'gaittypes');

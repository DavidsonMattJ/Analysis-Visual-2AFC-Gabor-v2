%j1_findpeaksinHeadPos.m

% updated for GABOR VERSION
% Here we will identify the peaks in head position data, for later 
% splitting  walk trajectories into individual gait cycles.

% Gabor v2 (UTS) VERSION


%% parameters:
% useful for debugging
% slightly slower, but prints a trial x trial summary of 
% head position (Y), with troughs overlayed.
visualizeResults=1;     

%add paths if new session/ running in isolation.
if ~contains(path, homedir)
set_myOnlineDirectories_vGaborv2;
end

cd(procdatadir)

% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);


Fs=60; % approx sampling rate - confirm with diff(Timevec) (was 90!)
%
% Note the threshold between peaks for detection might need to change per
% participant, based on e.g. step length.

%the int usetype in adjustGaitresultsBinned_2speed.m (called below)

%defaults:
pkdist=10; %samples between
pkheight=.005;


pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)

%%
for ippant =1:length(pfols)


    cd(procdatadir)
    
    %%load data from import job.
    load(pfols(ippant).name, 'HeadPos', 'trial_summaryTable', 'subjID', 'pnum');
    savename = pfols(ippant).name;
    disp(['Preparing j1 ' subjID]);


    % try to reorient to fig dir, make subfolder if absent
    try cd([figdir  filesep 'trial_headYpeaks' filesep 'p' pnum '_' subjID])
        pfigdir= pwd;
    catch
        mkdir([figdir  filesep 'trial_headYpeaks'  filesep 'p'  pnum '_' subjID]);
        cd([figdir  filesep 'trial_headYpeaks'  filesep 'p' pnum '_' subjID])
        pfigdir=pwd;
    end

    %% visualize the peaks and troughs in the head (Y) position.
    %Head position on Y axis, smooth to remove erroneous peaks.
    if visualizeResults
        figure(1); clf;
        set(gcf, 'units', 'normalized', 'position', [0.01,0.01, .9, .9], 'color', 'w', 'visible', 'off');
        %
        pcount=1; % plot indexer
        figcount=1; %figure counter.
    end



    %% some trials will/have been identified as bad, worth skipping in the below for-loop.
    
    badtrials = rejTrials_detectvGaborv2(subjID);

    for  itrial=1:size(HeadPos,1)
       
        %%
        if ismember(itrial,badtrials)
            disp(['Skipping bad trial ' num2str(itrial) ' for ' subjID]);
             
            %% if badtrial detected:            
                if visualizeResults
                    subplot(5,3,pcount)
                    plot(Timevec, trialD);
                    title('badtrial');

                    pcount=pcount+1;

                    if pcount==15

                        %print that figure, and reset trial count.
                        cd(pfigdir)
                        print('-dpng', [plotD.subjID ' trialpeaks ' num2str(figcount)]);
                        newfig=1;
                        figcount=figcount+1;
                        clf;
                        pcount=1;
                    end
                
                end
                %% continue after plotting.
                continue
        end % ismember bad trial
            
        
        
        % extract trial data:
        trialtmp = HeadPos(itrial).Y;
        trialD_sm = smooth(trialtmp, 5); % reduced smoothing factor (was 5).
        trialD = squeeze(trialtmp);


        Timevec = HeadPos(itrial).time;

        reptrial = find(trial_summaryTable.trial== itrial); % representative trialss

        if isempty(reptrial) % some trials may have missing data, if no targets were shown.
            continue
        end
        wlkSpeed = trial_summaryTable.walkSpeed(reptrial(1)); % 1 for slow, 2 for normal.
        
        if wlkSpeed==1
            wlkPrint = 'slow';
            minmaxlength_samps =[.15 .75].*Fs;

        else
            wlkPrint = 'normal';
            %flag gait durations outside this range (in sec) for review.
            minmaxlength_samps =[.25 .75].*Fs;

        end
        
        %we dont have targ in the frame data for AUD, instead, we can use
        %the trial summary info.


        %locations of peaks and troughs:
        [locs_ptr, locs_trtr]=deal([]);

        % extract click onsets for plots (remove omissions and Nans).
        trialRTs = trial_summaryTable.reactionTime(reptrial);
        %some have nan (due to detection in j0_import.
        trialClicks = trial_summaryTable.clickOnset(reptrial);
        clickRTs= trialRTs + trialClicks;
        clickRTs= clickRTs(~isnan(clickRTs));
        % FAs are basically the clicks recorded as 'too early'. now a Nan
        % in the reactionTime column.
        trialFAindexs = find(isnan(trial_summaryTable.reactionTime(reptrial)));
        if trialFAindexs
            % collect that click onset, based on index of nan:
            trialFAs= trial_summaryTable.clickOnset(reptrial(trialFAindexs));
            
        else
            trialFAs=nan;
        end


        % Note that we can skip when ppant was practicing/ stationary
        if  itrial>=5 % first X(4)? were stationary 

            %% start of peak / trough detection region.

            % some particular combinations of subject/trial,
            % need adjusted step lengths (hard coded after visual
            % inspection)

            % adjustGait_subjlist_vAUD;
          
            % if not a bad trial from visual inspection, continue.

            %find local peaks and troughs.

            [~, locs_p]= findpeaks(trialD_sm, 'MinPeakDistance',pkdist,'MinPeakProminence', pkheight);
            [~, locs_tr]= findpeaks(-trialD_sm, 'MinPeakDistance',pkdist ,'MinPeakProminence', pkheight);
            
            %for some participants, smoothing obscures thepeak, so revert
            %to original: (visual inspection required).
            % [~, locs_p]= findpeaks(trialD, 'MinPeakDistance',pkdist,'MinPeakProminence', pkheight);
            % [~, locs_tr]= findpeaks(-trialD, 'MinPeakDistance',pkdist ,'MinPeakProminence', pkheight);


            % extract nearest peaks and troughs from unsmoothed data:
            [~, locs_p_r] =  findpeaks(trialD, 'MinPeakDistance', pkdist, 'MinPeakProminence', pkheight);
            %find nearest in raw data, to peaks detected in smoothed version
            locs_ptr=zeros(1,length(locs_p));
            for ip=1:length(locs_p)
                [~, idx] = min(abs(locs_p_r - locs_p(ip)));
                locs_ptr(ip) = locs_p_r(idx);
            end

            %make sure no duplicates:
            locs_ptr= unique(locs_ptr);

            %same for troughs:
            [~, locs_tr_r] =  findpeaks(-trialD, 'MinPeakDistance', pkdist, 'MinPeakProminence', pkheight);

            %find nearest in raw data, to peaks detected in smoothed version
            locs_trtr=zeros(1,length(locs_tr));
            for ip=1:length(locs_tr)
                [~,idx] = min(abs(locs_tr_r - locs_tr(ip)));
                locs_trtr(ip) = locs_tr_r(idx);
            end

            %make sure no duplicates:
            locs_trtr= unique(locs_trtr);

            % for debugging: may want to plot the first result
                        % figure(10); clf;
                        % plot(Timevec, trialD,'b-');hold on;
                        % plot(Timevec, trialD_sm, 'k-');hold on;
                        % plot(Timevec(locs_ptr), trialD(locs_ptr), 'ro');
                        % plot(Timevec(locs_trtr), trialD(locs_trtr), 'bo');

            % end region for detection

            

            %% start region for clean/tidy peak deteciton
            % for stability, we want each trial to start and end in a trough.
            if locs_trtr(1) > locs_ptr(1) % if trial starts with peak.
                
                %remove that early peak.
                locs_ptr(1)=[];
                
                %insert trough, using minimum before first peak
                %                 [~, ftr ] = min(trialD(1:locs_ptr(1)));
                %                 locs_trtr = [ ftr,  locs_trtr];

            end
            %if trial ends with peak, add trough after at local minimum
            if locs_ptr(end) > locs_trtr(end) % if trial ends with peak.
                %insert trough, using minimum after last peak
                [~, ftr ] = min(trialD(locs_ptr(end):end));
                locs_trtr = [locs_trtr, locs_ptr(end)-1+ftr];
            end



            % finally, make sure that troughs and peaks alternate, if not, use
            % the maximum (peak) or minimum (trough) option
            % Should only be an issue at trial onset.
            % test if first peak is after a second trough , if so remove the
            % latter.

            %store copy to not mess with for loop indexing.
            newpks = locs_ptr;

            %%
            if length(locs_trtr) ~= length(locs_ptr)+1
                %% correct as necessary:

                % find doubled pks (most likely):
                for igait = 1:length(locs_trtr)-1
                    gstart = locs_trtr(igait);
                    gend = locs_trtr(igait+1);
                    try % inspect gait cycles for double peaks.
                        if locs_ptr(igait+1) < gend % if two peaks before gait cycle ends.
                            %retain the max height as true peak.
                            h1= trialD(locs_ptr(igait));
                            h2= trialD(locs_ptr(igait+1));
                            if h1>h2
                                %remove next peak
                                locs_ptr(igait+1)=[];
                            else %remove first peak.
                                locs_ptr(igait)=[];
                            end

                        end
                    catch
                    end

                end



                % now find doubled troughs

                for igait = 1:length(locs_ptr)
                    %for each peak, check only one trough before and after.
                    gstart = locs_trtr(igait);
                    try  gend = locs_trtr(igait+1);
                    catch
                        continue
                    end
                    if locs_ptr(igait) > gend % if two troughs before a peak.
                        %retain the min height as true trough.
                        h1= trialD(gstart);
                        h2= trialD(gend);
                        if h1>h2
                            %remove first trough
                            %                        deltroughs = [deltroughs, igait];
                            locs_trtr(igait)=[];
                        else %remove second trough.
                            locs_trtr(igait+1)=[];
                        end

                    end
                    %% check if this is the last peak, that there is only one trough remaining:
                    if igait == length(locs_ptr)
                        %last trough should be max, else error.
                        if locs_trtr(igait+1) ~= (locs_trtr(end))
                            % then retain only the minimum of the remaining
                            % troughs.
                            h1= trialD(locs_trtr(igait+1));
                            h2= trialD(locs_trtr(end));
                            if h1>h2
                                %remove first trough
                                %                        deltroughs = [deltroughs, igait];
                                locs_trtr(igait+1)=[];
                            else %remove second trough.
                                locs_trtr(end)=[];
                            end
                        end
                    end

                end % for each gait


            end % if more troughs than peaks.

            %% some sanity checks:
            % any very small/ big gaits?
            gaitszs = diff(locs_trtr);
            if any(gaitszs<minmaxlength_samps(1) ) || any(gaitszs>minmaxlength_samps(2) )
                %remove very short gaits:
                sml = find(diff(locs_trtr) <minmaxlength_samps(1)  );
                trscopy = locs_trtr;
                pkscopy=locs_ptr;
                while sml
                    %remove last peak and trough.
                    trscopy(sml+1)=nan;
                    pkscopy(sml)=nan;
                    %check remaining.
                    sml= find(diff(trscopy) <minmaxlength_samps(1)  );
                    disp('removed one small step');
                end
                locs_trtr= trscopy(~isnan(trscopy));
                locs_ptr= pkscopy(~isnan(pkscopy));

            end
            % make sure we go from trough to trough.
            if length(locs_trtr) ~= length(locs_ptr)+1

                disp(['!check gait count  trial:' num2str(itrial)]);

            end

            % any very small/ big gaits?
            gaitszs = diff(locs_trtr);
            if any(gaitszs>minmaxlength_samps(2) ) || any(gaitszs<minmaxlength_samps(1) )

                disp(['!check gait size trial ' num2str(itrial) ': ' subjID]);

            end
            %%


            %add these peaks and troughs to trial structure data.
            HeadPos(itrial).Y_gait_peaks = locs_ptr;
            HeadPos(itrial).Y_gait_peaks_sec = Timevec(locs_ptr);
            HeadPos(itrial).Y_gait_troughs = locs_trtr;
            HeadPos(itrial).Y_gait_troughs_sec = Timevec(locs_trtr);




        end
        %
        relvrows = find(trial_summaryTable.trial==itrial);
        trialTarg = trial_summaryTable.targOnset(relvrows);
        trialRTs= trial_summaryTable.clickOnset(relvrows);
        % for visualisation, plot the targ and rt onsets.

        %we may want to vis the trial onsets and response:
        if visualizeResults

            % small function to plot per panel:
            plotD=[];
            plotD.itrial=itrial;
            plotD.figcount=figcount;
            plotD.hy= trialD;
            plotD.trg = trialTarg;
            plotD.clk= trialRTs;
            plotD.Timevec= Timevec;
            plotD.troughs = locs_trtr;
            plotD.pks = locs_ptr;
            plotD.pfigdir= pfigdir;
            plotD.walkPrint= wlkPrint;
            plotD.subjID= subjID;
            plotD.FAs= trialFAs;
            % see funtion defn below.
            newFig=quickplot(pcount,plotD);
            if newFig==1
                pcount=2;
                figcount=figcount+1;
            else
                pcount=pcount+1;

            end

        end

    end %itrial.
    %% after all trials, print the final figure (in case uneven subplots/trial counts).
    if visualizeResults
        cd(pfigdir)
        print('-dpng', [subjID ' trialpeaks ' num2str(figcount)]);
    end
    %resave with new data in structure.
    cd(procdatadir)
    save(savename, 'HeadPos', '-append');

end % ippant

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calling this function above.
function newfig= quickplot(pcount, plotD)
% function to handle subplot placement, and printing/counting new figures.
%% newfig is a flag to start the next figure, toggles subplot count.
figure(1);
newfig=0;
if plotD.itrial==16 || pcount ==16
    %print that figure, and reset trial count.

    cd(plotD.pfigdir)
    print('-dpng', [plotD.subjID ' trialpeaks ' num2str(plotD.figcount)]);
    newfig=1;
    clf;
    pcount=1;
end
subplot(5,3,pcount);
plot(plotD.Timevec, plotD.hy   );
hold on;
%add extracted peaks and troughs.
plot(plotD.Timevec(plotD.pks), plotD.hy(plotD.pks), 'or');
plot(plotD.Timevec(plotD.troughs), plotD.hy(plotD.troughs), 'ob')
ylabel('Head position');
xlabel('Time (s)');
title(['Trial ' num2str(plotD.itrial) '-' plotD.walkPrint]);
yrange = range(plotD.hy);
ylim([(min(plotD.hy) - yrange) max(plotD.hy)]);

yyaxis right
for itarg= 1:length(plotD.trg);
    hold on;
    plot([plotD.trg(itarg), plotD.trg(itarg)],[0, 1], 'k-');
end
% plot clickRTs?
for iclk = 1:length(plotD.clk)
    % plot(plotD.Timevec, plotD.clk, 'r');
    plot([plotD.clk(iclk) plotD.clk(iclk)], [0 1], 'r-');
end

%plot FAs
for iFAs = 1:length(plotD.FAs)
    plot([plotD.FAs(iFAs) plotD.FAs(iFAs)], [0 1], 'm-');
end

ylabel('targ-click');
ylim([0 2]);

% axis tight
end


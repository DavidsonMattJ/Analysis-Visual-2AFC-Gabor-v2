% j2C_append_doubGCData_gaborv2
% building off the back of j2- > adding an extra column for the
% classification of an event in either a LRL or RLR step sequence.
% necessary for showing double step proportions of DVs etc.

% append event percentages, relative to stride.

%previously was extracting for both LRL and RLR progression, now including
%just sequential steps in stride.
% e.g.  step 1-2-3 is stride 1; step 3-4-5 is stride two, etc.
% ensures no double dipping of onset events.

%laptop:
%add paths if new session/ running in isolation.
if ~contains(path, homedir)
    set_myOnlineDirectories_vGaborv2;
end

cd(procdatadir)
% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%
%%


for ippant= 1:nsubs
    cd(procdatadir)    %%load data from import job.
    load(pfols(ippant).name, ...
        'HeadPos', 'trial_summaryTable', 'subjID');
    savename = pfols(ippant).name;
    disp(['Preparing j2 cycle data... ' savename]);

    %% Gait extraction.
    % per event, find the relevant gait. then add the gPcnt relative
    % to LRL anr RLR sequence (if possible).

    allevents = size(trial_summaryTable,1);

    count0_1_100=zeros(2,3,1); % debug counter for edge cases, increments within ppant.

    for ievent= 1:allevents
        %guardclause to skip trial if  practice or stationary, or flagged
        %for rejection.

        itrial = trial_summaryTable.trial(ievent);

        badtrials=rejTrials_detectvGaborv2(subjID); %toggles skip based on bad trial ID
        skip=0;
        if ismember(itrial,badtrials)
            skip=1;
        end

        if trial_summaryTable.walkSpeed(ievent)==0 ||skip
            %>>>>> add this new info to our table,
            %trgOnset:
            trial_summaryTable.respO_gPcnt_strideinTrial(ievent)= nan;
            trial_summaryTable.trgO_gPcnt_strideinTrial(ievent)= nan;

            continue
        else


            %for each event, convert the target onset time, and reaction onset
            %time, to a (%) of a double gait cycle (i.e. percent in
            %stride).

            tTrial = trial_summaryTable.trial(ievent);

            trs = HeadPos(tTrial).Y_gait_troughs;
            trs_sec = HeadPos(tTrial).Y_gait_troughs_sec;
            pks = HeadPos(tTrial).Y_gait_peaks;
            trialTime = HeadPos(tTrial).time;
            eventcolumns = {'targOnset', 'clickOnset'};
            savecols = {'trgO', 'respO'};

            for ieventtype=1:2
                %%Target onsets first:
                % find the appropriate gait. Final trough that the event is later
                % than:
                thisEvent = trial_summaryTable.([ eventcolumns{ieventtype} ])(ievent);

                gindx = find(thisEvent>trs_sec,1,'last');


                %Determine if we should allocate to the first or second
                %step in a stride, based on trial progression.
                %%
                allsteps= trs_sec; % timestamp of all troughs.

                % to ensure no L/R order effects, start the stride allocation on
                % alternate steps per trial.

                if mod(tTrial,2)==0 % even numbers
                    %shift by one.
                    allsteps= allsteps(2:end); % omitting the first trough.
                end

                %set up an array that shows which stride we are on.
                [stepAllocation,strideAllocation]=deal(nan(1,length(allsteps)));
                stridecounter=1;
                stepinStride=1; % will be 1 or 2, alternating.
                stepsfrom= repmat([1,2], [1,100]);

                for istep= 1:length(allsteps)

                    strideAllocation(istep)= stridecounter; % which stride #in trial?
                    stepAllocation(istep)= stepsfrom(istep); % first or second step in stride?

                    %after every three heel strikes, increment stride counter
                    if mod(istep,2)==0 % even numbers
                        stridecounter=stridecounter+1;
                    end

                end
                %debug:
                % disp([allsteps, stepAllocation', strideAllocation'])

                %%

                %don't fill for misses, or before/after first/last step in a trial.
                if thisEvent ==0 || isempty(gindx) || gindx== length(trs)  %
                    
                    trial_summaryTable.([savecols{ieventtype}  '_gPcnt_strideinTrial'])(ievent)= nan;


                else
                    %%  extract info we need:

                    %Base the reference based on trial progression.
                    %which stride in trial, and which step in that stride?

                    % gindx indexes into trs_sec (all troughs). For even
                    % trials, allsteps dropped the first trough, so the
                    % allocation arrays are shifted by 1 relative to gindx.
                    if mod(tTrial, 2) == 0
                        gindx_alloc = gindx - 1;
                    else
                        gindx_alloc = gindx;
                    end

                    % Guard: event falls before the first trough included
                    % in this trial's stride scheme (only possible on even
                    % trials when gindx == 1).
                    if gindx_alloc < 1 || gindx_alloc > length(stepAllocation)
                        trial_summaryTable.([savecols{ieventtype} '_gPcnt_strideinTrial'])(ievent) = nan;
                        continue
                    end

                    tmpStrideAllocation = strideAllocation(gindx_alloc);
                    tmpstepinStride = stepAllocation(gindx_alloc); % was it after the first or second heelstrike?

                    if tmpstepinStride==1 % event was within first step of stride.
                        try trialstrideSamps = trs(gindx):trs(gindx+2)-1;
                        catch
                            trialstrideSamps=nan;
                            disp('Warning! event occuring outside of trial data: ');
                        end
                    elseif tmpstepinStride==2 % event was within second step of stride.

                        try trialstrideSamps = trs(gindx-1):trs(gindx+1)-1;
                        catch
                            trialstrideSamps=nan;
                            disp('Warning! event occuring outside of stride data:')
                        end

                    end



                    gaitSamps = trialstrideSamps;
                    doubgPcnts=[]; % store output.

                    if  any(isnan(gaitSamps))
                        continue % pass over.
                    end

                    gaitTimes = HeadPos(tTrial).time(gaitSamps);

                    %  take as proportion of total.
                    tDur = gaitTimes(end)- gaitTimes(1);
                    tE= thisEvent-gaitTimes(1);
                    gPcnt= round((tE/tDur)*100);

                    % note that the result can now be 0 ! which
                    % messes with later calculations.
                    if ismember(gPcnt, [0,1,100])
                        idx=([0,1,100]==gPcnt);
                        count0_1_100(ieventtype, idx,1)=count0_1_100(ieventtype,idx,1)+1 ;

                        % flip coin and give to either 1 or 100.
                        % if randi([0 1])
                        % gPcnt=1;
                        % else
                        % gPcnt=100;
                        % end
                    end


                    doubgPcnts= gPcnt;


                    %>>>>> add this new info to our table,                    
                    trial_summaryTable.([savecols{ieventtype}  '_gPcnt_strideinTrial'])(ievent)= doubgPcnts;



                end % if not misses
            end % targ and response onset


            %%       % plot to be safe:
            % if ismember(itrial, [20,21,22,23,24,25])
            % 
            % 
            %     clf; subplot(1,3,1:2)
            %     plot(HeadPos(tTrial).time, HeadPos(tTrial).Y);
            %     title(num2str(tTrial))
            %     %add troughs
            %     hold on;
            %     plot(HeadPos(tTrial).time(trs), HeadPos(tTrial).Y(trs), ['ob'], 'LineWidth',2 )
            %     shg
            %     %add targets
            %     trialrows = find(trial_summaryTable.trial==tTrial);
            %     tOns = trial_summaryTable.targOnset(trialrows);
            %     yyaxis right
            %     for itrg= 1:length(tOns)
            % 
            %         plot([tOns tOns], [0 mean(HeadPos(tTrial).Y)], 'r-','LineWidth',2)
            % 
            %     end
            %     ylim([0 5]);
            %     % now shade the strides
            %     stridecols= repmat({'b', 'r'}, [1,20]);
            %     for istep=1:length(trs_sec)-1;
            %         xvec = [trs_sec(istep), trs_sec(istep) trs_sec(istep+1) trs_sec(istep+1)];
            % 
            %         yvec=[0 .25 .25 0];
            % 
            %         ph=patch(xvec, yvec, 'k', 'FaceAlpha', 0.2/stepAllocation(istep));
            %         text(mean(trs_sec(istep:istep+1)), .25, num2str(stepAllocation(istep)), 'VerticalAlignment','top', 'fontsize', 10)
            % 
            %         % and add the strides.
            %         yvec = [0.25 .5 .5 0.25];
            % 
            %         ph=patch(xvec, yvec, stridecols{strideAllocation(istep)}, 'FaceAlpha', 0.2);
            %         text(mean(trs_sec(istep:istep+1)), .5, num2str(strideAllocation(istep)), 'VerticalAlignment','top', 'fontsize', 10)
            % 
            %     end
            %     text(0, .25, 'step')
            %     text(0, .5, 'stride')
            %     % show the final distirbution of onsets by trial stride:
            %     subplot(1,3,3)
            %     histogram(trial_summaryTable.trgO_gPcnt_strideinTrial(trialrows), length(trialrows));
            % 
            %     shg
            %     shg
            %     % showt
            % end

        end % if not practice etc


    end % each row in table (event)

    %% quick sanity check:
    %
    if visualiseOnsets
    allgO = trial_summaryTable.trgO_gPcnt;
    allgO = allgO(~isnan(allgO)); % remove nans
    allgO_doubinTrial = trial_summaryTable.trgO_gPcnt_strideinTrial; % 
    allgO_doubinTrial= allgO_doubinTrial(~isnan(allgO_doubinTrial )); % removenans
    
    %store the counts per g pcnt, plot both. 
    mycounts= nan(2,100);
    for ig=1:100
        mycounts(1,ig) = length(find(allgO==ig));
        mycounts(2,ig) = length(find(allgO_doubinTrial==ig));
    end
    % plot both.
    clf
    subplot(211);
    bar(1:100, mycounts(1,:));
    title(['target counts for single step (all), N=' num2str(sum(mycounts(1,:)))]);
    subplot(212);
    hold on;
    bar(1:100, mycounts(2,:));
    title(['target counts for tride (two steps - all), N=' num2str(sum(mycounts(2,:)))])

    shg
    end


    
    %%
    disp(['Finished appending DOUBLE gait percentage data for ... ' subjID]);
    save(savename, 'trial_summaryTable','-append');
end % participant


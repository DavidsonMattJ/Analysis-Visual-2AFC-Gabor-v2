% j2_appendgaitcycleData_wSDT_vGaborv2

% Here we will load the summary data table, and add event markers  new columns, for when events
% happened relative to % gait cycle completion.
% For now, events of interest are target onset, and response onset. 

% jobs:
% -  append event percentages
% -  appends whether L/R ft (now in j2a).


% Gabor discrim version2 (UTS).

visualiseOnsets = 1; % print participant level onsets (histogram) for sanity check.

%add paths if new session/ running in isolation.
if ~contains(path, homedir)
set_myOnlineDirectories_vGaborv2;
end

cd(procdatadir)


%% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr);



%%

for ippant= 1:nsubs
    cd(procdatadir)    %%load data from import job.
    load(pfols(ippant).name, ...
        'HeadPos', 'trial_summaryTable', 'subjID');
    savename = pfols(ippant).name;
    disp(['Preparing j2 cycle data... ' savename]);

    %% Gait extraction.
    % Per trial (and event), extract gait samples (trough to trough), normalize along x
    % axis, and store various metrics.

    count0_1_100=zeros(2,3,1); % debug counter for edge cases,
    % will increment.
    % dimsare (targ/resp, [0,1,100], count).

    allevents = size(trial_summaryTable,1);



    for ievent= 1:allevents
        %guardclause to skip trial if  practice or stationary
        %trials of experiment.

        %or if previously identified as a bad trial :
            itrial = trial_summaryTable.trial(ievent);               
            badtrials=rejTrials_detectvGaborv2(subjID); %toggles skip based on bad trial ID
            skip=0;
            if ismember(itrial,badtrials)
                skip=1;
            end

        if trial_summaryTable.walkSpeed(ievent) ==0 || skip
            %>>>>> add this new info to our table,
            %trgOnset:
            trial_summaryTable.trgO_gCount(ievent) = nan;
            trial_summaryTable.trgO_gPcnt(ievent)= nan;
            trial_summaryTable.trgO_gDur(ievent)= nan;
             trial_summaryTable.trgO_gStart_sec(ievent)= nan;
              trial_summaryTable.trgO_gFin_sec(ievent)= nan;
               trial_summaryTable.trgO_gStart_samp(ievent)= nan;
                trial_summaryTable.trgO_gFin_samp(ievent)= nan;
               

            
            %response
            trial_summaryTable.respO_gCount(ievent) = nan;
            trial_summaryTable.respO_gPcnt(ievent)= nan;
            trial_summaryTable.respO_gDur(ievent)= nan;
               trial_summaryTable.respO_gStart_sec(ievent)= nan;
              trial_summaryTable.respO_gFin_sec(ievent)= nan;
               trial_summaryTable.respO_gStart_samp(ievent)= nan;
                trial_summaryTable.respO_gFin_samp(ievent)= nan;
      
              

            continue
        else
            

            %For each event, conver the target onset time, and reaction onset
            %time, to a % of a gait cycle.

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

                %don't fill for misses, or before/after first/last step in a trial.
                if thisEvent ==0 || isempty(gindx) || gindx== length(trs)  %
                    trial_summaryTable.([savecols{ieventtype}  '_gCount'])(ievent) =nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gPcnt'])(ievent)= nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gDur'])(ievent)= nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gStart_sec'])(ievent) = nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gFin_sec'])(ievent) = nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gStart_samp'])(ievent) = nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gFin_samp']) (ievent)= nan;

                else % extract info we need:
                
                    
                    % extract the event as % gait.[0 100]
                    gaitsamps =trs(gindx):trs(gindx+1)-1; % avoid counting edge bins twice. (-1?)
                    
                    gaitTimes = HeadPos(tTrial).time(gaitsamps);
                    
                    %! [old way], but stretching may have introduced noise
                    %at boundaries.
                    % resizeT= imresize(gaitTimes', [1,100]);
                    % gPcnt = dsearchn(resizeT', thisEvent);
                    
                    %! [new way]
                    %  take as proportion of total. 
                    tDur = gaitTimes(end)- gaitTimes(1);
                    tE= thisEvent-gaitTimes(1);
                    gPcnt= round((tE/tDur)*100);
                    
                   %% keep track of edge cases for later debugging:                                      
                   if ismember(gPcnt,[0,1,100])
                       %store a count for later debugging.
                       idx=([0,1,100]==gPcnt);
                       count0_1_100(ieventtype, idx,1)=count0_1_100(ieventtype,idx,1)+1 ;
                       
                       % % flip coin and give to either 1 or 100.
                       % warning, introduces noise. This isn't necessary if
                       % we avoid counting edge bins twice, as above in 
                       % gaitsamps ...

                       % if randi([0 1])
                       %     gPcnt=1;
                       % else
                       %     gPcnt=100;
                       % end
                   end



                            %% debug plot.                            
                           %  clf; 
                           %  plot(HeadPos(tTrial).time, HeadPos(tTrial).Y);
                           %  hold on;
                           %  plot(HeadPos(tTrial).time(trs), HeadPos(tTrial).Y(trs), 'bo');
                           %  %and targets:
                           %  %all tO
                           %  yyaxis right
                           %  tindex = find(trial_summaryTable.trial==tTrial);
                           %  ttO = trial_summaryTable.targOnset(tindex);
                           % plot(ttO', ones(1, length(ttO)), 'ro','linew',2)
                           % 
                           % trO = trial_summaryTable.clickOnset(tindex); %trial response
                           % plot(trO', ones(1, length(trO)), 'bx', 'linew',2)
                           % shg
%%

                    %>>>>> add this new info to our table,
                    trial_summaryTable.([savecols{ieventtype}  '_gCount'])(ievent) = gindx; %trgO_gcount, or respO_gcount.
                    trial_summaryTable.([savecols{ieventtype}  '_gPcnt'])(ievent)= gPcnt;
                    trial_summaryTable.([savecols{ieventtype}  '_gDur'])(ievent)= round(gaitTimes(end)-gaitTimes(1),3);
                    trial_summaryTable.([savecols{ieventtype}  '_gStart_sec'])(ievent) = [gaitTimes(1)];
                    trial_summaryTable.([savecols{ieventtype}  '_gFin_sec'])(ievent) = [gaitTimes(end)];
                    trial_summaryTable.([savecols{ieventtype}  '_gStart_samp'])(ievent) = [gaitsamps(1)];
                    trial_summaryTable.([savecols{ieventtype}  '_gFin_samp'])(ievent) = [gaitsamps(end)];

                    if gaitsamps(1)<1
                        disp('debug')
                    end

                end % if not misses
            end % targ and response onset

        end % if not practice

    end % each row in table (event)

% debug. Show trgOs
% ptrgOs = trial_summaryTable.trgO_gPcnt;
% ptrgOs = ptrgOs(~isnan(ptrgOs));
% hist(ptrgOs, 100);

if any(trial_summaryTable.respO_gPcnt ==0)|| any(trial_summaryTable.trgO_gPcnt ==0)
% check code
error('0 gait index for event has been retained, check code');
end
    
    %% include zscored version of RTs:
    allRTs = trial_summaryTable.reactionTime;
    userow = find(~isnan(allRTs));
    zRTs = zscore(allRTs(userow));
    %place in table as new columnL
    trial_summaryTable.z_clickRT = nan(size(trial_summaryTable,1),1);
    trial_summaryTable.z_clickRT(userow) = zRTs;



   %% critical for later stages, remove the data for 'bad' trials.
   badtrials= rejTrials_detectvGaborv2(subjID); % pulled from rejTrials_detectvGabor.
   allts = trial_summaryTable.trial;
   remtrials = ismember(allts, badtrials);
   trial_summaryTable(remtrials,:)=[];

%% quick sanity check.
if visualiseOnsets
%%
%all
    allgO = trial_summaryTable.trgO_gPcnt;
    %also slow and norm (only)
    slowts = find(trial_summaryTable.walkSpeed==1);
    normts = find(trial_summaryTable.walkSpeed==2);
    slowgO= trial_summaryTable.trgO_gPcnt(slowts);
    normgO= trial_summaryTable.trgO_gPcnt(normts);

%remove nans (else shown at indx 0 on histogram).
allgO= allgO(~isnan(allgO));
slowgO= slowgO(~isnan(slowgO));
normgO= normgO(~isnan(normgO));

clf;
pData={slowgO, normgO, allgO};
titlesare={'slow', 'natural','combined'};
for id=1:3
    tmpD= pData{id};
    subplot(2,2,id); % top row histograms.
    %histogram seems buggy at edge cases.
    mycounts= nan(1,100);
    for ig=1:100
        mycounts(ig) = length(find(tmpD==ig));
    end
    bar(1:100, mycounts);
title(titlesare{id})
    
end
subplot(2,2,4) % below plot bar_count
bar(count0_1_100(1,:)); grid on;
set(gca,'XTickLabel', {'0','1','100'});
xlabel('total targ counts at 0,1,100')
    shg
end
%%
   disp(['Finished appending gait percentage data for ... ' subjID]);
   save(savename, 'trial_summaryTable','-append');
end % participant


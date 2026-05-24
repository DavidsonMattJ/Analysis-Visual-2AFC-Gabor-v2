% Discrimination experiment (Gabor)-
%%  Import from csv. FramebyFrame, then summary data.

%%%%%% v4: Gabor data (Discrimination version) %%%%%%


cd(datadir)

%list all participant files:

pfolsframe = dir([pwd filesep '*framebyframe.csv']);

pfolssummary= dir([pwd filesep '*_trialsummary.csv']);

% show ppantfiles/ numbers in command window
tr= table((1:length(pfolssummary))',{pfolssummary(:).name}' );
disp(tr)


sancheckplots=0; % toggle to show sanity checks for data flags (FA per trial etc).


%% Per csv file, import and wrangle into Matlab Structures, and data matrices:


for ippant =1:length(pfolsframe)
    cd(datadir)

    %% load subject data as table.
    filename = pfolsframe(ippant).name;
    %extract name&date from filename:
    ftmp = find(filename =='_');
    subjID = filename(1:ftmp(1)-1);
    %%
    
    pnum = sprintf('%02d', ippant);
    savename= ['p' pnum '_' subjID  '_summary_data'];


    %read table
    cd(datadir);
    opts = detectImportOptions(filename,'NumHeaderLines',0);
    disp('reading large frame x frame file now...');
    T = readtable(filename,opts);
    % ppant = T.participant{1};
    disp(['Preparing participant ' pnum]);

  
    %% First create HeadPos structure used in later scripts:

    alltrials = unique(T.trialNumber); % note the index starts at 0, so shift.
           
    ntrials= length(alltrials);
    %preallocate a struct array with each field initialised to  []
    HeadPos= struct('X', cell(ntrials,1),...
        'Y', cell(ntrials,1),...
        'Z', cell(ntrials,1),...
        'time', cell(ntrials,1),...
        'blockType', cell(ntrials,1));
    
    % duplicate structure for EyePos 
    EyePos=HeadPos;
%%

    for itrial= 1:ntrials
        
        rowIndx = find(T.trialNumber==alltrials(itrial));
        %store:
        HeadPos(itrial).Y =  T.head_Y(rowIndx);
        HeadPos(itrial).X =  T.head_X(rowIndx);
        HeadPos(itrial).Z =  T.head_Z(rowIndx);
        HeadPos(itrial).time =  T.trialTime(rowIndx);
        
        
        %we also have eye pos
        EyePos(itrial).gazeDir_X= T.gazeDirection_X(rowIndx);
        EyePos(itrial).gazeDir_Y= T.gazeDirection_Y(rowIndx);
        EyePos(itrial).gazeDir_Z= T.gazeDirection_Z(rowIndx);
        
        % copy over the gaze data.
        EyePos(itrial).gazeHit_X= T.gazeHit_X(rowIndx);
        EyePos(itrial).gazeHit_Y= T.gazeHit_Y(rowIndx);
        EyePos(itrial).gazeHit_Z= T.gazeHit_Z(rowIndx);
        
        EyePos(itrial).gazeOrigin_X= T.gazeOrigin_X(rowIndx);
        EyePos(itrial).gazeOrigin_Y= T.gazeOrigin_Y(rowIndx);
        EyePos(itrial).gazeOrigin_Z= T.gazeOrigin_X(rowIndx);

        EyePos(itrial).gazeHitObject = T.gazeHitObject(rowIndx);
        EyePos(itrial).gazeAngularSpeed = T.gazeAngularSpeed(rowIndx);
        EyePos(itrial).time =  T.trialTime(rowIndx);

    end


    %% save before next step:
    cd(procdatadir)


    try save(savename, 'HeadPos', 'EyePos', 'subjID', 'pnum', '-append');
    catch
        save(savename, 'HeadPos' ,'EyePos','subjID', 'pnum');
    end
    %
    %% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ now summary data



    cd(datadir)
    
    nsubs= length(pfolssummary);
    %
    filename = pfolssummary(ippant).name;

    %extract name&date from filename:
    ftmp = find(filename =='_');
    subjID = filename(1:ftmp(end)-1);
    %read table
    opts = detectImportOptions(filename,'NumHeaderLines',0);
    
    T = readtable(filename,opts);
    rawSummary_table = T;

    disp(['Preparing participant (summary) ' pnum]);


    
    %% extract Target onsets per trial (as table).
    % and Targ RTs, other stats from the csv if relevant.
    frametrials = length(alltrials);
    alltrials = unique(T.trial);
    
  
    
    disp([subjID  ' has ' num2str(length(alltrials)) '(summary)'])


    if frametrials~=length(alltrials)
        disp([ subjID ': trial mismatch!'])
        % when this is the case, it is likely because of the 'extra' trials
        % bug in the frame x frame data. This can now be removed.

        maxTrialdata = length(alltrials);
        HeadPos= HeadPos(1:maxTrialdata);
    end

    %we want to repair the table as we go.
    trial_summaryTable= T;

    %remove some cols we don't need:
    trial_summaryTable.date=[];

    try trial_summaryTable.qStep=[]; % not all have it for some reason.
    catch
    end

    %rename correctResponse to targCor, and tonePresent to signalPresent  (match this/later scripts).
    colid1= find(contains(trial_summaryTable.Properties.VariableNames, 'correctResponse'));
    trial_summaryTable.Properties.VariableNames{colid1}= 'targCor';
    
    %change strings to binary for signal present/absent.

    sigP = find(contains(trial_summaryTable.signalPresent,'True'));
    sigA = find(contains(trial_summaryTable.signalPresent,'False'));
    %wipe first.
    trial_summaryTable.signalPresent=[];
    trial_summaryTable.signalPresent(sigP)=1;
    trial_summaryTable.signalPresent(sigA)=0;


    %% There's a bug.
    % Some trials start with a RT on frame 1 (first 11ms). It's due to
    % requiring both triggers to be pressed to start a trial.
    % Some participants don't release straight away.
    % Remove these data rows.

    a= find(trial_summaryTable.reactionTime < .1);
    b= find(trial_summaryTable.reactionTime>0); % no responses are -10. (so this finds those RTs betwen 0 and .1
    c= find(trial_summaryTable.trialID==0); % restrict to only those at trial onset.
    tmp = intersect(a,b);

    remrows = intersect(tmp,c);
    %remove those rows.
    trial_summaryTable(remrows,:)=[];
    
    %% note the remaining list are worth addressing, as RTs are faster than
    %plausible (< .1s)
    queryrows = setdiff(tmp,c);

    

    %% now work through all the data per trial (a single walk),
    % and commit the data to a new structure.
    for itrial= 1:length(alltrials)
        
        relvrows = find(trial_summaryTable.trial==alltrials(itrial)); % unity idx at zero.

        %create vectors for storage:
        tOnsets = trial_summaryTable.targOnset(relvrows);
        tRTs = trial_summaryTable.reactionTime(relvrows);
        tCor = trial_summaryTable.targCor(relvrows);


        %% Reaction time region (tidy / reclassify some values).

        
        %note that negative RTs, indicate that no response was recorded:
        tOmit = find(tRTs<=0);
        if ~isempty(tOmit)
            %             tCor(tOmit) = NaN; % don't count those incorrects, as a mis identification.
            tRTs(tOmit)=NaN; % remove no respnse
        end

        
        %% we also want to reject very short RTs (reclassify as a FA).

        if any(find(tRTs<0.15))
            % disp(['Trial: ' num2str(thistrial) ' Suspicious RT']);
            checkRT= find(tRTs<.15);
            %debug to check:

            % reclassify data (mark as incorrect, and record as FA).
            for ifalseRT = 1:length(checkRT)
                indx = relvrows(checkRT(ifalseRT));
                
                trial_summaryTable.reactionTime(indx)=nan;
                trial_summaryTable.correctResponse(indx) =  nan;
                
            end
            %
            % update for next catch:
            tOnsets = trial_summaryTable.targOnset(relvrows);
            tRTs = trial_summaryTable.reactionTime(relvrows);

        end
        %%
        % seems some FA are missing, do a quick check to reclassify double
        % responses and determine FAs.
        % frequently multiple triggers are pulled (could be self
        % corrections?).
        % find the *second* trigger within a window, and convert to a FA.
        % these are identfiable as duplicate entries to the same target onset.
        if any(diff(tOnsets)==0)


            if sancheckplots
                %plot
                clf;
                %head data
                plot(trialInfo(itrial).times, HeadPos(itrial).Y); hold on;;

                yyaxis right
                %trg and rts:
                for itrgo = 1:length(tOnsets)
                    plot([tOnsets(itrgo) tOnsets(itrgo)], [0, 1], 'k-o')
                    text(tOnsets(itrgo), 1, num2str(itrgo))
                end
                for irt = 1:length(tRTs)
                    plot([tRTs(irt) tRTs(irt)], [0, 1], 'r-o')
                    text(tRTs(irt), 1, num2str(irt))
                end

                ylim([0 2]);
                shg
            end




            % reclassify data (mark as incorrect, and record as FA).
            checkRT= find(diff(tOnsets)==0);
            %debug to check:

            % reclassify data (mark as incorrect, and record as FA).
            for ifalseRT = 1:length(checkRT)
                indx = relvrows(checkRT(ifalseRT));
                %repair table to avoid counting twice!
                % note we are assuming the second response is the error:
                repindx= indx+1;

                trial_summaryTable.FA_rt(repindx) =  trial_summaryTable.targRT(repindx);
                trial_summaryTable.clickRT(repindx) = NaN;
                trial_summaryTable.targOnset(repindx)=NaN;
                trial_summaryTable.signalPresent(repindx)= NaN;
                trial_summaryTable(repindx, 11:16)= table(NaN);
            end
            % update for next catch:
            tOnsets = trial_summaryTable.targOnset(relvrows);
            tRTs = trial_summaryTable.targRT(relvrows);
        end  %any duplicates

        % end RT region


    end
    %repair trial indexing:
    trial_summaryTable.trial =  trial_summaryTable.trial +1;
    trial_summaryTable.block =  trial_summaryTable.block +1;
    trial_summaryTable.trialID =  trial_summaryTable.trialID +1;

    
    
    % change all -10 to nan (these were no response recorded.
    chRT= find(trial_summaryTable.reactionTime==-10);
    trial_summaryTable.reactionTime(chRT)=nan;

    % add a column for easy categorisation (H, M, FA, CR).= 1,2,3,4
    %
    presData = trial_summaryTable.signalPresent;
    sigPresent = find(trial_summaryTable.signalPresent==1);
    sigAbsent= find(trial_summaryTable.signalPresent==0);
    respPres = find(trial_summaryTable.targResponse ==1);
    respAbs = find(trial_summaryTable.targResponse ==0);

    Hits = intersect(sigPresent, respPres);
    Misses = intersect(sigPresent, respAbs);
    FalseAlarms = intersect(sigAbsent, respPres);
    CorrectRejections = intersect(sigAbsent, respAbs);

    trial_summaryTable.SDTcat(Hits) = 1;
    trial_summaryTable.SDTcat(Misses) = 2;
    trial_summaryTable.SDTcat(FalseAlarms) = 3;
    trial_summaryTable.SDTcat(CorrectRejections) = 4;

    %save this structure for later analysis per gait-cycle:
    disp(['Saving trial summary data ... ' subjID]);
    cd(procdatadir)

    save(savename, 'trial_summaryTable','-append');


end % participant
% end % with out without.

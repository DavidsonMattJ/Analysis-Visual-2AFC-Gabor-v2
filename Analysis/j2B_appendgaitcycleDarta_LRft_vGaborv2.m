% j2B_appendgaitcycleData_LRft_vGaborv2

% Here we will load the summary data table, and add event markers  new columns, for when events
% happened relative to % gait cycle completion.
% For now, events of interest are target onset, and response onset.

% jobs:
% -  appends whether L/R ft (now in j2B).


% % % % Gabor v2 task VERSION UTS% % %

%laptop:
%add paths if new session/ running in isolation.
if ~contains(path, homedir)
    set_myOnlineDirectories_vGaborv2;
end
cd(procdatadir)
%% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%
%%


for ippant =1:nsubs
    cd(procdatadir)    %%load data from import job.
    load(pfols(ippant).name, ...
        'HeadPos',  'trial_summaryTable', 'subjID');
    savename = pfols(ippant).name;
    disp(['Preparing j2 cycle data... ' savename]);

    %% Gait extraction.
    % Per trial (and event), extract gait samples (trough to trough), normalize along x
    % axis, and store various metrics.

    allevents = size(trial_summaryTable,1);



    for ievent= 1:allevents
        %guardclause to skip trial if  practice, stationary or flagged for rejection.

        itrial = trial_summaryTable.trial(ievent);

        badtrials=rejTrials_detectvGaborv2(subjID); %toggles skip based on bad trial ID
        skip=0;
        if ismember(itrial,badtrials)
            skip=1;
        end


        if trial_summaryTable.walkSpeed(ievent)==0 || skip

            %>>>>> add this new info to our table,
            trial_summaryTable.trgO_gFoot(ievent) = {'nan'};
            trial_summaryTable.respO_gFoot(ievent) = {'nan'};

            continue
        else



            %for each event, determine whether the gait cycle is Left-Right
            %ft, or Right-Left ft. This is done by observing the sway (left
            % to right), in head position.


            % Head position data:
            tmpPos=  squeeze(HeadPos(itrial).Y);
            tmpSway = squeeze(HeadPos(itrial).Z);
            tmpwalkDir = squeeze(HeadPos(itrial).X);

            %is walk direction (x axis) increasing or decreasing?
            if mean(diff(tmpwalkDir))<0
                walkDir= 'decreasing';
            elseif mean(diff(tmpwalkDir))>0
                walkDir= 'increasing';
            else
                error('debug!');
            end

            trialTime = HeadPos(itrial).time;

            %% sanity  check:
            % plotD={tmpPos, tmpSway};
            % ylabs={'Y-height', 'Z-sway'};
            % for idata=1:2
            %     figure(1)
            %     subplot(2,2,1+ 1*(idata-1));
            %     plot(tmpwalkDir, plotD{idata});
            %     xlabel('raw X pos');
            %     ylabel(ylabs{idata})
            %     title(['trial ' num2str(itrial) ', X is ' walkDir])
            %     % note that matlab will always rearrange the x-axis to be neg to positive!
            %     %so if decreasing, adjust:
            %     if strcmp(walkDir, 'decreasing')
            %         set(gca,'XDir', 'reverse');
            %     end
            % 
            % 
            %     subplot(2,2,3 +1*(idata-1));
            %     plot(trialTime,  plotD{idata});
            %     xlabel('trial time');
            %     ylabel(ylabs{idata})
            % end
            % shg
            %%
            % quick classification:
            if strcmp(walkDir, 'increasing')
                % Then more positive z values
                % are left side of the body.

                Zpos = 'LHS'; % allocentric
                Zneg= 'RHS';

            elseif strcmp(walkDir, 'decreasing')
                % the return trajectory), pos z values are RHS
                Zpos = 'RHS';
                Zneg='LHS';

            end


            % note that for each event (row int the table), we could have
            % separate gaits for target onset and response.
            % So calculate both
            savecols = {'trgO', 'respO'};
            for ieventtype=1:2

                gStart = trial_summaryTable.([savecols{ieventtype} '_gStart_samp'])(ievent);
                gEnd = trial_summaryTable.([savecols{ieventtype} '_gFin_samp'])(ievent);

                % if no response though, continue (this is the case for
                % respO when missing - no response = nan.
                if any(isnan([gStart gEnd]))
                    disp('debug');

                    %update table: 
                    trial_summaryTable.([savecols{ieventtype} '_gFoot'])(ievent) = {'nan'};
                    continue
                else
                    % is the z value increasing or decreasing, relative to feet placement?

                    %
                    midlineS = mean(tmpSway([gStart, gEnd])); % start and end (return) of sway.
                    meanS = mean(tmpSway(gStart:gEnd));

                    gaitSway = tmpSway(gStart:gEnd);

                    if strcmp(walkDir, 'increasing') && meanS>midlineS
                        % Zpos= LHS, Zneg = RHS
                        ft='LR'; % Zpos (vs midline),  .:. leaning left
                    
                    elseif strcmp(walkDir, 'increasing') && meanS<midlineS
                        % Zpos= LHS, Zneg = RHS

                        ft= 'RL'; %Zneg (vs midline), .:. leaning right

                    elseif strcmp(walkDir, 'decreasing') && meanS>midlineS
                        % Zpos= RHS, Zneg = LHS                        
                        
                        ft= 'RL'; %Zpos (vs midline) .:. leaning right
                    elseif strcmp(walkDir, 'decreasing') && meanS<midlineS
                        % Zpos= RHS, Zneg = LHS                        
                        ft= 'LR';%Zneg (vs midline) .:. leaning left
                    end


                    %% sanity check:
                    % figure(2)
                    % clf;
                    % 
                    % subplot(311);
                    % plot(trialTime, tmpSway); %
                    % xlabel('Time');
                    % title(['Trial ' num2str(itrial) ' raw sway - Z axis']);
                    % subplot(312);
                    % plot(tmpwalkDir, tmpSway); %
                    % if strcmp(walkDir, 'decreasing')
                    %     set(gca,'XDir', 'reverse');
                    % end
                    % 
                    % subplot(313)
                    % plot(trialTime, tmpSway, 'color', [.8 .8 .8]);
                    % 
                    % % % to ease interp, reorient to allocentric:
                    % % % problem is +- Z values are L/R side
                    % % % of body, based on walking direction.
                    % %
                    % % if strcmp(walkDir, 'OUT')
                    % %     set(gca,'Ydir','reverse')
                    % %      ylabel('Right -Left')
                    % % else
                    % %     ylabel(ylab)
                    % % end
                    % 
                    % %overlay
                    % hold on;
                    % plot(trialTime(gStart:gEnd), gaitSway, 'k', 'linew', 2)
                    % hold on;
                    % plot(trialTime([gStart, gEnd]), tmpSway([gStart, gEnd]), 'ro')
                    % text(trialTime(gStart), tmpSway(gStart), ft(1), 'color', 'b', 'fontsize', 15)
                    % text(trialTime(gEnd), tmpSway(gEnd), ft(2), 'color', 'b', 'fontsize', 15)
                    % %add text to be sure
                    % title(['Walking ' walkDir ', gait: ' ft])
                    % shg

                    %%
                    %>>>>> add this new info to our table,


                    trial_summaryTable.([savecols{ieventtype} '_gFoot'])(ievent) = {ft};
                end
            end % eventtype (targetonset or resp in gait).

        end % if not practice

    end % each row in table (event)

    disp(['Finished appending gait percentage data for ... ' subjID]);
    save(savename, 'trial_summaryTable','-append');
end % participant


% plotPFX_forInspection_vGaborv2.m
%
% Per-participant inspection plots.
%
% Layout (2 rows × 2 cols, one PNG per participant):
%   staircase, acc, HR, FAR
%
% Data source: trial_summaryTable
%
% Prerequisites:
% job.plot_ppantResults:                J0_ImportVRData
% job.plot_distribution_outliers:
%
% Saves PNGs to: figdir/Participant inspection for exclusion/

close all;
% -------------------------------------------------------------------------
%  Directories (guard: skip if homedir already on path)
% -------------------------------------------------------------------------
if ~contains(path, homedir)
    set_myOnlineDirectories_vGaborv2;
end

% =========================================================================
%  TOGGLE FLAGS
% =========================================================================
jobs=[];
jobs.plot_ppantResults=1;
jobs.plot_distribution_outliers=0;

% -------------------------------------------------------------------------
%  File-name suffix — define locally so script can run standalone
% -------------------------------------------------------------------------
if ~exist('usenewTrialStride', 'var')
    usenewTrialStride = 1;   % default: use '_new' stride-epoch files
end
if usenewTrialStride
    appendp = '_new';
else
    appendp = '';
end

%% Prep/Load data ------------------------------------------------------------------
cd(procdatadir );

pfols =dir([pwd filesep 'p_*']);
nppants= length(pfols);



%% Plot configuration ---------------------------------------------------------
pCols       = {'k', 'b', 'm'};       % one colour per DV
fontsize    = 15;
speedLabels = {'Slow', 'Normal', 'Combined'};
speedCols= {'b', 'r', 'b','r'};
DVnames = {'Acc', 'HR', 'FAR'};     % DVs to display (edit to add/remove)
ncols   = length(DVnames) + 1;            % DV columns + rightmost Fourier R² column

% Configure the figure window once — clf (called per participant below)
% clears axes but preserves these figure-level properties across iterations.
figure(1);clf;
set(gcf, 'Color', 'w', 'Units', 'normalized', 'Position', [0 0 .5 1]);

%% Output directory -----------------------------------------------------------
inspectDir = fullfile(figdir, 'Participant inspection for exclusion');
if ~isfolder(inspectDir)
    mkdir(inspectDir);
end
if jobs.plot_ppantResults
%% Per-participant loop --------------------------------------------------------
for ippant = 1:nppants

    clf   % clear axes for this participant; figure window properties are preserved

    %% Pre-collect data for all speeds and DVs -----------------------
   cd(procdatadir)
   %DV data is in grandAvg, but staircase is in trialsummary table

   load(pfols(ippant).name, 'trial_summaryTable', 'subjID', 'pnum');
   psubjID=subjID;

   %staircase plot will be accuracy over time. 
   slowts = find(trial_summaryTable.walkSpeed==1);   
   normts = find(trial_summaryTable.walkSpeed==2);
   allts = sort([slowts;normts]); % both (excludes stationary/practice).

   % signal present only subsets  
   sigP = find(trial_summaryTable.signalPresent==1);
   slowts_hr = intersect(slowts, sigP);
   normts_hr = intersect(normts, sigP);

   sigA = find(trial_summaryTable.signalPresent==0); % used below.

   % plot both staircases overlaid (accum accuracy - should converge to
   % threshold set in params).
   
   %can calc grandaverage Acc, HR, FAR now:
   speedlists = {slowts, normts, allts};
   plotDVs= nan(3,3); % DVs x 3, speeds x 3
   for ispeed = 1:3 % slow, matural, combined
    trialsubset = speedlists{ispeed};

    %acc = 
    thisAcc = mean(trial_summaryTable.targCor(trialsubset), 'omitnan');
    %HR calculated on only signal present trials: 
    sigPtrials = intersect(sigP, trialsubset);
    thisHR = mean(trial_summaryTable.targCor(sigPtrials), 'omitnan');

    %FAR calculated on only signal absent trials: 
    sigAtrials = intersect(sigA, trialsubset);
    thisFAR = mean(trial_summaryTable.targResponse(sigAtrials), 'omitnan'); % 

    %now store: 
    plotDVs(1,ispeed) = thisAcc;
    plotDVs(2,ispeed) = thisHR;
    plotDVs(3,ispeed) = thisFAR;
   end

%% progress to plotting, staircases first:

   subplot(2,2,1);
   staircols={'b','r','b','r'};
   stairlns={'-','-',':',':'};
   stairtrials= {slowts, normts, slowts_hr, normts_hr};
   
   for idata= 1:4
       
       tdata = trial_summaryTable.targCor(stairtrials{idata});
       %remove nans (non responses.
       tdata = tdata(~isnan(tdata));
       t_cumsum = cumsum(tdata)';
        
       t_rolAv = t_cumsum./[1:length(t_cumsum)];
       hold on;
       plot(t_rolAv, 'color', staircols{idata}, 'LineStyle', stairlns{idata});

       % add final. 
       text(length(t_cumsum), t_rolAv(end),  sprintf('%.2f',t_rolAv(end)),...
           'HorizontalAlignment','left', ...
           'VerticalAlignment','bottom', ...
           'Color', speedCols{idata}, ...
           'Fontsize', fontsize*1.5);


   end
   legend('slow','natural', 'slow (hr)', 'natural (hr)', 'autoupdate','off');
   
   title('Staircase')
   ylim([0.45, 1]); 
   ylabel('cumulative average');
   xlabel('Target count');
   set(gca,'fontsize', fontsize)
   hold on; yline(0.5, 'k:','LineWidth',2)
   
   shg
   
    %% Draw subplots for DVs -----------------------------------------------------------
    for idv = 1:length(DVnames)

            %% Bar chart: DV at each speed --------------------
            subplot(2,2, 1+idv)

            pDV = plotDVs(idv,:);
            bar(1:3, pDV, 'FaceColor', pCols{idv}, 'FaceAlpha', .5);
            hold on;
            
            for ib=1:3
            text(ib, pDV(ib), sprintf('%.2f',pDV(ib)),...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','bottom','FontSize', fontsize)
            end

            ylabel(DVnames{idv})
            
            title( [psubjID '(' num2str(ippant) '), ' DVnames{idv} ], 'FontSize', fontsize);
            xlabel('Speed','FontSize', fontsize);
            set(gca, 'FontSize', fontsize, 'Xticklabel', {'slow','normal','combined'});
            ylim([0 1]);
shg
    end % idv
    shg

    

    %% Figure-level title -----------------------------------------------------
    sgtitle(['Participant ' num2str(ippant) ': ' psubjID], ...
        'FontSize', fontsize, 'Interpreter', 'none');

    %% Save PNG ---------------------------------------------------------------
    print(fullfile(inspectDir, ...
        ['Participant ' pnum ' ' psubjID ' Basics per speed' appendp]), ...
        '-dpng');
    shg

end % ippant
end % job. 

if jobs.plot_distribution_outliers
%% here plot the distribution of our main DVs, to determine outliers.

% wrangle
%% for each, plot the histogram and id outliers.
figure(2);clf;
set(gcf, 'Color', 'w', 'Units', 'normalized', 'Position', [0 0 .5 1]);
speedCols={'b','r','k'};

 
plotDVs=[];
for ippant=1:nppants
    for idv=1:length(DVnames)
        plotDVs(idv,ippant,1) = GFX_grandAvg(ippant,1).(['grand_' DVnames{idv}]); % Acc, HR, FAR, (3 speeds, slow, normal, combined).
        plotDVs(idv,ippant,2) = GFX_grandAvg(ippant,2).(['grand_' DVnames{idv}]); % Acc, HR, FAR, (3 speeds, slow, normal, combined).
        plotDVs(idv,ippant,3) = GFX_grandAvg(ippant,3).(['grand_' DVnames{idv}]); % Acc, HR, FAR, (3 speeds, slow, normal, combined).
    end
end

for idv = 1%:3%length(DVnames)

            % histogram per.
            for ispeed=1:3
            
            distributionDV = plotDVs(idv,:,ispeed);
            outls= find(isoutlier(distributionDV));
            subplot(length(DVnames),3, ispeed + length(DVnames)*(idv-1))
            histogram(distributionDV,40);
            

            xlabel(DVnames{idv})
            ylabel('count')
            
            
            % identify outliers in legend?
            ylim([0 5])
            xlim([0 1])
            if ~isempty(outls)
                %build outlier string.
                outstring=[];
                for io= 1:length(outls)
                outstring= [outstring ' '  subjIDs{outls(io)}];

                end

                title({[speedLabels{ispeed} ', ' DVnames{idv} ];...
                    ['outliers are: ' outstring]});
            
            else % normal title.
                title([speedLabels{ispeed} ', ' DVnames{idv} ]);
            end
            % ylim([0 1]);
            end

% print output for debugging: 
[srtDV,srtDVindx]=sort(distributionDV);

%% show in command window
tr=table(srtDVindx',srtDV',{subjIDs{srtDVindx}}' );
tr.Properties.VariableNames{1}= 'ippant';
tr.Properties.VariableNames{2}= [speedLabels{ispeed} '_' DVnames{idv} ];
tr.Properties.VariableNames{3}= 'subjID';

disp(tr)

end % idv
    %%
end
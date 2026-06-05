% function job_MSplot_grandaverageSummary(GFX)
% this is a (close) to MS ready script, focusing on differences in
% av step duration, Acc, and RT only.

%panel 1: avStepDur
%panel 2: Acc
%panel 3: RT

omitbadPpants=1; 

cd(procdatadir); cd('GFX');
load('GFX_grandAvg_data.mat',...
    'GFX_grandAvg', 'GFX_headY',...
    'avStepDuration');

speedCols  = [0 0 1; ...
    1, .7 ,.25];


figure(1);
set(gcf,'Units', 'normalized', 'Position', [.6 .4 .4 .4 ], 'color', 'w')
set(gcf,'Units', 'normalized', 'Position', [0 0 1 1 ], 'color', 'w')
clf

fsize=20;

% no omit to start with.
% omitPpants=[22,6,15,5,28,2, 4,16,     17,20];
omitPpants=[]; % exclusions applied..

plotGFX= GFX_grandAvg;
plotStepDurs = avStepDuration;

if omitbadPpants        
    plotGFX(omitPpants,:)= [];
    plotStepDurs(omitPpants,:)=[];

end

nppants= size(plotGFX,1);
% First Adding top panels for step and stride time-series.
gcfields= {'gc', 'doubgc'};
cyclenames = {'step', 'stride'};
for igc=1:2
    for ispeed= 1:2

        if igc==1
            subplot(2,3,1)
        else
            subplot(2,3,2:3)
        end
        hold on
        %
        plotHead=[];
        for ippant=1:nsubs
        plotHead(ippant,:) = GFX_headY(ippant,ispeed).([gcfields{igc} '_resamp']);
        end
    
        %add shaded: 
        stE= CousineauSEM(plotHead);
        gM= mean(plotHead,1,'omitnan');

        sh=shadedErrorBar(1:length(gM), gM, stE, {'color', speedCols(ispeed,:)});

        xlabel(['Head position as % of ' cyclenames{igc} '-cycle']);
        ylabel('detrended head height (m)')
        
        
        formatGaitAxis(gca,1:length(gM))
        set(gca,'fontsize',fsize)
box on
        
        
    end
end
shg
%%


% first wrangle the data for plots: 
% dvfields={'acc', 'HR', 'FAR', 'RT','dprime','crit'};
dvfields={'Acc', 'rt'};
dvfields={'HR', 'FAR'};
dvfields={'dprime', 'crit'};

allboxData={}; % step, acc, rt

allboxData{1}= plotStepDurs;
for iDV= 1:2

%wrangle the data.

slowD=  [plotGFX(:,1).(['grand_' dvfields{iDV}])];
natD=  [plotGFX(:,2).(['grand_' dvfields{iDV}])];
tmpdata= [slowD', natD'];

if strcmp(dvfields{1}, 'Acc')
    tmpdata= tmpdata.*100;
end
allboxData{iDV+1}= tmpdata;
end

ylabels={'Step Duration (sec)', 'Accuracy (%)', 'Reaction time (sec)'};
ylimsare= [0 1;...
    0 100;...
    0 1];

% now plot each:
for iplot = 1:3
%% plots box, whiskers, and scatter, 
% data = points x groups
% 
% params.scatSize;
% params.cols
% params.scatCol
% params.scatAlpha=.1
% params.boxlinewidth
% params.plotScatter
% params.boxtype= 'boxplot'

% cla
subplot(2,3,iplot+3)
params=[];
params.scatSize=40;
params.cols  = [0 0 1; ...
    1, .7 ,.25];
% params.scatCol='b';
params.boxlinewidth=2;
params.scatAlpha= .2;
params.plotScatter=1;
params.boxtype= 'boxplot'; % or 'boxplot'
params.boxwidth=.8;
params.joinScatter=1;
h=box_and_scatter(allboxData{iplot},params);
% xlim([.8 2.2])
% adjust some aparms.
% ylim(ylimsare(iplot,:));
my_scaleYrange(.9, allboxData{iplot})
hold on;
% yline(.5,'--','color','k', 'LineWidth', 2)
%
ylabel(ylabels{iplot});
title([ylabels{iplot}]);
xlabel('Walking speed');
set(gca,'XTickLabels', {'slow', 'natural'},'fontsize', fsize);


%add stats:
ttestdata = allboxData{iplot};
[h,pval,ci,stat] = ttest(ttestdata (:,1), ttestdata (:,2));
ppos = 'bottom'; % vertalign.
if pval <.001
    pstat= '***';
elseif pval < .01
    pstat = '**';
elseif pval < .05
    pstat = '*';
else
    pstat = 'ns';
    ppos = 'bottom'; %veralign.
end
hold on;
%add bar between (at 80% plot height) and signifiy significance level:
ysat = get(gca,'ylim');
pplace = ysat(1) + (ysat(2)- ysat(1)).*.9;
plot([1,2], [pplace,pplace], 'k-');
text(1.5, pplace, pstat, 'VerticalAlignment',ppos, 'HorizontalAlignment','center', 'fontsize',fsize*2);



shg
end % iplot
%% print('-dpng', ['Group Calibration summary N(' num2str(nppants) ') '  groupname{igroup} ]);
shg
%%

% end


function formatGaitAxis(gca, xvec)
% formatGaitAxis  Apply standard gait-cycle axis tick formatting.
%   Labels start (0%), midpoint (50%), and end (100%) of the gait cycle.
midp = xvec(ceil(length(xvec) / 2));
set(gca, 'FontSize', 15, ...
    'XTick',      [xvec(1), midp, xvec(end)], ...
    'XTickLabel', {'0%', '50%', '100%'});
end
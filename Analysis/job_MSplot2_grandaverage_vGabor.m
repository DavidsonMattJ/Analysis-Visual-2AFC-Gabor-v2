% function job_MSplot2_grandaverage(GFX)
%% quick plot, inspect grand effects (not resolved over stride)

%panel 1: Accuracy
%panel 2: HR
%panel 3: FAR
%panel 3: RT
%panel 4: Dprime
%panel 5: criterion

omitbadPpants=0; 

cd(procdatadir); cd('GFX');
load('GFX_grandAvg_data.mat',...
    'GFX_grandAvg', 'GFX_headY',...
    'avStepDuration');
%%
% dvfields={'acc', 'HR', 'FAr', 'RT','dprime','crit'};
dvfields={'Acc', 'HR', 'FAR', 'rt','dprime','crit'};

ymax = [1,1,1,1,2.5,1.5];

clf

fsize=15;

% no omit to start with.
% omitPpants=[22,6,15,5,28,2, 4,16,     17,20];
omitPpants=[];

plotGFX= GFX_grandAvg;
% plotStepDurs = avStepDuration;

if omitbadPpants        
    plotGFX(omitPpants,:)= [];
    plotStepDurs(omitPpants,:)=[];

end

nppants= size(plotGFX,1);

for iDV= 1:6

%wrangle the data.


slowD=  [plotGFX(:,1).(['grand_' dvfields{iDV}])];
natD=  [plotGFX(:,2).(['grand_' dvfields{iDV}])];
tmpdata= [slowD', natD'];


subplot(2,3,iDV);
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
params=[];
params.scatSize=40;
params.cols  = [0 0 1];
% params.scatCol='b';
params.boxlinewidth=3;
params.scatAlpha= .2;
params.plotScatter=1;
params.boxtype= 'boxplot'; % or 'boxplot'
params.boxwidth= .4;

h=box_and_scatter(tmpdata,params);

% adjust some aparms.
ylim([0 ymax(iDV)])
%
ylabel(dvfields{iDV});
title([dvfields{iDV} ' n=' num2str(nppants)]);
set(gca,'XTickLabels', {'slow', 'natural'},'fontsize', fsize);


%add stats:

[h,pval,ci,stat] = ttest(tmpdata(:,1), tmpdata(:,2));
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
pplace = ysat(1) + (ysat(2)- ysat(1)).*.8;
plot([1,2], [pplace,pplace], 'k-');
text(1.5, pplace, pstat, 'VerticalAlignment',ppos, 'HorizontalAlignment','center', 'fontsize',fsize);



shg
end % iDV
%% print('-dpng', ['Group Calibration summary N(' num2str(nppants) ') '  groupname{igroup} ]);
shg
%%

% end
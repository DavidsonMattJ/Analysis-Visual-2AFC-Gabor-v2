% function job_MSplot_Headtimeseries(GFX)
% this is a (close) to MS ready script, focusing on differences in
% head height during slow and natural pace, showing the resampling in
% effect

set_myOnlineDirectories_vGaborv2;

cd(procdatadir);
cd('GFX');
load('GFX_grandAvg_data.mat',...
    'GFX_headY');

speedCols  = [0 0 1; ...
    1, .7 ,.25];

%defaults changed from j1_ due to the reduced range of time series:
%WAS
% pkdist=10; %samples between
% pkheight=.005;
%NOW
pkdist=5; %samples between
pkheight=.00001;


figure(1);
set(gcf,'Units', 'normalized', 'Position', [.6 .4 .4 .4 ], 'color', 'w')
set(gcf,'Units', 'normalized', 'Position', [0 0 1 1 ], 'color', 'w')
clf

fsize=20;

% no omit to start with.

omitPpants=[]; % exclusions applied..

striderange=[400,50]; % min samps we expect the stride trough to be (slow,nat)
nppants= size(GFX_headY,1);
% First Adding  panels for step/stride time-series.
% per ppant, raw overlayed with low alpha.

gcfields= {'gc', 'doubgc'};
cyclenames = {'step', 'stride'};
datanames= {'raw', 'resamp'};
useMean=0; % 0 for overlay with shade, 1 for mean.
igc=2; 
for idata=1:2
    for ispeed= 1:2

        if idata==1
            subplot(2,2,1)
        else
            subplot(2,2,2)
        end
        hold on
        %
        plotHead=[];
        for ippant=1:nsubs
            plotHead(ippant,:) = GFX_headY(ippant,ispeed).([gcfields{igc} '_' datanames{idata}]);
        end

        if idata==2%useMean==1
            %add shaded:
            stE= CousineauSEM(plotHead);
            gM= mean(plotHead,1,'omitnan');

            sh=shadedErrorBar(1:length(gM), gM, stE, {'color', speedCols(ispeed,:)});
            xlabel([cyclenames{igc} '-cycle (%)']);
            ylabel('detrended head height (m)')

            formatGaitAxis(gca,1:length(gM))

        else
            % overlay.
            % per participant, plot with alpha, and stop at last reasonable
            %% point.
            % clf
            maxend=[];
            for ippant = 1:nsubs

                tmpTS = plotHead(ippant,:);
                %% sanity check, plot turning points (troughs).
                [~, locs_tr]= findpeaks(-(smooth(tmpTS,5)), 'MinPeakDistance',pkdist,'MinPeakProminence', pkheight);
                % % limit to
                % plot(tmpTS); hold on;
                % plot(locs_tr, tmpTS(locs_tr), 'ro')
                % shg50
                stoptrough = find(locs_tr>(striderange(igc)));
                % if 


                % end
                tmpTS_trim = tmpTS(1:locs_tr(stoptrough(1))); % go to first or second trough

                % ph=plot((1:length(tmpTS_trim))./60, tmpTS_trim, 'color', [speedCols(ispeed,:), .5]); % add alpha
                ph=plot((1:length(tmpTS))./60, tmpTS, 'color', [speedCols(ispeed,:), .5]); % add alpha
                hold on
                %%
                maxend(ippant)=locs_tr(stoptrough(1)); % for x axis
            end
            
            % xts= get(gca,'XTick');
            % set(gca,'Xtick')
            xlabel(['Time (sec)']);
            shg

        end % overlau
      
        set(gca,'fontsize',fsize)
        % box on


    end
end
shg
%%



function formatGaitAxis(gca, xvec)
% formatGaitAxis  Apply standard gait-cycle axis tick formatting.
%   Labels start (0%), midpoint (50%), and end (100%) of the gait cycle.
midp = xvec(ceil(length(xvec) / 2));
set(gca, 'FontSize', 15, ...
    'XTick',      [xvec(1), midp, xvec(end)], ...
    'XTickLabel', {'0%', '50%', '100%'});
end
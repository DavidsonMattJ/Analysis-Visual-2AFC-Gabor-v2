% set up stepwise selection 
% STEPWISE_FOURIER_SELECTION
% Finds the best subset of time-series for a fourier1 group-level fit
%
% INPUTS:
%   X          - [40 x T] matrix of time series (rows = series, cols = time)
%   xvec       - time vector, e.g. 1:10
%   min_series - minimum series to retain (default: 30)
%   hz_range   - [lo hi] acceptable Hz window (default: [1.5 2.5])
%   predropped   - indices to exclude before selection begins (default: [])
%                  e.g. [1 2] removes series 1 and 2 from all consideration
%

%
% OUTPUTS:
%   results    - struct with best subset, R² trajectory, Hz estimates, flags
%
% EXAMPLE:
%   results = stepwise_fourier_selection(X, 1:10, 30, [1.5 2.5], [1,2,3]);


normON=1;
normtype='relchange';
% HR for slow/natural stride
dataplot=[];
for ippant= 1:size(GFX_TargPosData,1);
dataplot(ippant,:) = GFX_TargPosData(ippant,2,3).doubgc_binned_HRCalc; %dims (ippant, speed, foot)
end
if normON
dataplot = applyNorm(dataplot,normtype);
end

pidx1  = ceil(linspace(1, 101, 11));  % single step: 10 bins
  mdiff = round(mean(diff(pidx1)) / 2);
  xvec  = pidx1(1:end-1) + mdiff;


  %% Gabor task (nat)
  results_r2= stepwise_fourier_selection(dataplot,xvec, 34, [1.8, 2.2], ...
[3] );  %suggests %5, 30, 27, 16
%%
  stepwise_fourier_selection_composite(dataplot,xvec,35, [1.8,2.2],[])
  %%  
  %% Gabor (nat)
  results_hz= stepwise_fourier_selection_hz(dataplot,xvec, 35, [1.9, 2.2], ...
[] ); %starts in range  
  
  %% Gabor task (slow)
  
   results_hz= stepwise_fourier_selection_hz(dataplot,xvec, 32, [1.8, 2.2], ...
[16,6,3] ); % % biggest drop for p16,6,4
   %% .:. increase R2
   results_hz= stepwise_fourier_selection(dataplot,xvec, 35, [1.8, 2.2], ...
[16,6,3] ); %starts in range!
  %%
  


function data = applyNorm(data, normtype)
% applyNorm  Per-participant normalisation of a [nSubs × nBins] matrix.
%
%   Each row is normalised relative to that participant's row mean, so
%   the normalised grand mean across bins is approximately zero (or one).
%
%   normtype options:
%     'absolute'  – subtract row mean
%     'relative'  – divide by row mean, minus 1  (percent change from mean)
%     'relchange' – (x – mean) / mean             (classic relative change)
%     'normchange'– (x – mean) / (x + mean)       (bounded symmetric change)
%     'db'        – 10·log10(x / mean)             (decibel change)
pM       = mean(data, 2, 'omitnan');
meanVals = repmat(pM, 1, size(data, 2));
switch normtype
    case 'absolute',   data = data - meanVals;
    case 'relative',   data = data ./ meanVals - 1;
    case 'relchange',  data = (data - meanVals) ./ meanVals;
    case 'normchange', data = (data - meanVals) ./ (data + meanVals);
    case 'db',         data = 10 * log10(data ./ meanVals);
    otherwise
        warning('applyNorm: unknown normtype ''%s'' — returning raw data.', normtype);
end
end
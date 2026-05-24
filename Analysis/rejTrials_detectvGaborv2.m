function badtrials = rejTrials_detectvGaborv2(subjID)

% certain bad trials (identified in figure folder).
% look at 'TrialHeadYPeaks', to see if the tracking dropped out on any
% trials. Or if participants were not walking smoothly.

% Abbreviations for rejected trials:
% s = poor signal quality (drop-outs/discontinuities)
% g= poor gait extraction (head tracking unclear).

% Step in to reject particular ppant+trial combos

badtrials=nan; % used as a flag in case subjID is incorrect. 

switch subjID

   case '10006' % p01
        badtrials = []; % none.


    case '10111' %p02; ! Revisit, poor gait extraction on early slow trials.
        
        badtrials=[22,68,93,100,201];
    
    
    case '10159' % p03
        badtrials= []; % none
    
    
    case '10207' %p04
        badtrials=[]; %none

    case '10270' %p05
        badtrials=[]; %none

    case '10297'
        badtrials= 102; %signal

    case '10381'
        badtrials=[]; %none
    
    case '10402'
        badtrials=[]; %none
    case '10441'
        badtrials=[]; %none
    case '10444'
        badtrials=[]; %none
    case '10465'
        badtrials=[]; %none
    
    case '10471'
        badtrials=[]; %none
    case '10477'
        badtrials=[25]; % gait
    case '10537'
        badtrials=[]; %none
    case '10543'
        badtrials=[9,31:33,... % signal
                 50,58,59];  % gait

    case '10555'
         badtrials = [22,81,... %gait
                     86,90]; %signal
        
    case '10561'
        badtrials=[]; %none

    case '10573'
         badtrials=103; %signal
    
    case '10642'
        badtrials= [4,5,... % missing?
            81]; %gait
    
    case '10669'
        badtrials= []; %none
    
    case '10672'
    badtrials= []; %none
    
    case '10837'  %p22
        badtrials = []; %none

    case '10846'
        badtrials = []; % none

    case '10855'
        badtrials = []; % none

    case '10891'
        badtrials= []; % none
    
    case '10897' % p26 %REJECT?
        badtrials =  []; % none, but slow gait is noisy
    
    case '9517'
        badtrials = []; %none
    
    case '9598'
        badtrials= []; %none
    
    case '9643'
        badtrials= 145; %signal
            
    case '9649'
        badtrials= 175; %signal
    case '9670'
        badtrials =[]; %none

    case '9736'
        badtrials= [5,21]; %gait
    case '9769'
        badtrials = []; %none
    case '9850'
        badtrials= [111,188,195]; %gait
    case '9871'
        badtrials = []; %none

    case '9904'
        badtrials= [130,131,134,136,195,196]; % gait
    case '9955'
        badtrials= []; % none
    case '9979'
        badtrials =[]; %none
    
    
    case '9982'  %p39  REJECT - shuffling slow gait.
        badtrials= [5,21:30,...] % gait... many more.
            14]; % signal
    
    case '9988' %p40 REJECT?
        badtrials= [61,62,135,... % gait
            126,140,141,176:178,212:220]; % signal

end % 

%%
if isnan(badtrials)
    disp('NO subjID')
end



%%
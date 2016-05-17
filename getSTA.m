function [sta,sem] = getSTA(stimMat,spbool,winSize,varargin)
%% function [sta,sem] = getSTA(stimMat,spbool,winSize,[scaleTGL],[spikeOrder],[excludeHigherOrder])
% ==============================================
% Inputs:
%       stimMat:
%           is an N x M matrix where N is the number of time points (each being 1
%           ms). Each column is a stimulus to compute the STA for.
%       spbool:
%           N x 1 vector of zeros or ones where each row is a timepoint
%           (correspoding to stimMat) and a one indicates that a spike
%           occured at that time.
%       winSize:
%           the number of milliseconds before and after the spike to
%           average
%       [scaleTGL]:
%           Binary to either perform feature scaling of the inputs or not.
%           Defaults to yes (1)
%       [spikeOrder]:
%           Changes the type of spike to find the STA for
%           1 = singlets
%           2 = doublets... Defaults to 1
%       [excludeHigher]:
%           binary to either exclude higher order events(1) or keep them
%           (0). e.g. excludeHigher = 1 would find 0 1 0 as a singlet, but 0 1 1 0 as
%           not two singlets. Defaults to 0
% ======================================
% Outputs:
%       sta:
%           a 2(winSize)+1 x M matrix of the spike triggered average for
%           all stimulus dimensions
%       sem:
%           a 2(winSize)+1 x M matrix of the SEM of the spike triggered average for
%           all stimulus dimensions
% ========================================
%% varargin handling
numVargs = length(varargin);
if numVargs > 3
    error('too many optional input arguments')
end

optargs = {1,1,0};
optargs(1:numVargs) = varargin;
[scaleTGL,spikeOrder,excludeHigher] = optargs{:};
%% scale if desired
if scaleTGL
    for ii = 1:size(stimMat,2)
        stimMat(:,ii) = (stimMat(:,ii)-nanmean(stimMat(:,ii)))./(nanstd(stimMat(:,ii)));
    end
end

%% find events
matchStr = ones(1,spikeOrder);
if excludeHigher
    matchStr = [0 matchStr 0];
end
times = strfind(spbool',matchStr);
if isempty(times)
    fprintf('No events of order %i found\n',spikeOrder)
    sta = [];
    sem = [];
    return
end


%% Find triggered values
triggered = zeros(winSize * 2 + 1, size(stimMat,2),length(times));

for ii = 1:length(times) %loop over every spike
    
    % if a spike is found st the window overlaps the trial start, pad
    % the beginning with nans and get the triggered stimulus
    if times(ii)-winSize<1
        temp = nan(winSize * 2 + 1,size(stimMat,2));
        temp(end-times(ii)-winSize +1:end,:) = stimMat(1:(times(ii) + winSize),:);
        triggered(:,:,ii) = temp;
        % If a spike is found st the window overlaps the trial end, pad the
        % end with nans and get the triggered stimulus
        
    elseif (times(ii)+winSize)>size(stimMat,1)
        temp = nan(winSize * 2 + 1,size(stimMat,2));
        temp(1:winSize + size(stimMat,1)-times(ii) + 1,:) = stimMat(times(ii) - winSize:end,:);
        triggered(:,:,ii) = temp;
        % Get the triggered stimulus.
        
    else
        triggered(:,:,ii) = stimMat((times(ii) - winSize):(times(ii) + winSize),:);
    end
    
    
end

%% compute mean and S.E.M. across all spikes
sta = nanmean(triggered,3);
sem = nanstd(triggered,0,3)./sqrt(length(times));

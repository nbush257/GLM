function [sta,se] = getSTA(spikes,Xin,winSize,varargin)
%% function sta = getSTA(matrix,winSize,[spikeTypeFlag])
% ==============================================
% Inputs:
%       matrix:
%           is an N x 4 matrix where N is the number of time points (each being 1
%           ms). The first column is the binary spike train, the remaining
%           columns are FX,FY,M respectively
%     
%       winSize:
%           the number of milliseconds before and after the spike to
%           average
%       [spikeTypeFlag]:
%           Changes the type of spike to find the STA for.
%           's' = singlets
%           'seh' = singlets exclude higher order
%           'd' = doublets
%           'deh' = doublets exclude...
%           't' = triplets
%           'teh' = triplets exclude...
%           'q' = quadruplets
%           'qeh' = quadruplets exclude...
% ======================================
% Outputs:
%       trigFX: the STA of FX (the second column of the input matrix)
%       trigFY: the STA of FY (the third column of the input matrix)
%       trigM: the STA of M (the fourth column of the input matrix)
% ========================================
% Error bars are standard error of the mean

spikeTypeFlag = '';
if length(varargin)==1
    spikeTypeFlag = varargin{1};
end

if ~ismember(spikeTypeFlag,{'s','seh','d','deh','t','teh','q','qeh'})
    spikeTypeFlag = 's';
    disp('automatically using singlets')
end
if iscolumn(spikes)
    spikes = spikes';
end

s = strfind(spikes,1);
seh = strfind(spikes,[0 1 0]);seh = seh +1;
d =strfind(spikes,[1 1]);
deh = strfind(spikes,[0 1 1 0]);deh = deh+1;
t = strfind(spikes,[1 1 1]);
teh = strfind(spikes,[0 1 1 1 0]);teh = teh+1;
q = strfind(spikes,[1 1 1 1]);
qeh = strfind(spikes,[0 1 1 1 1 0]);qeh = qeh+1;

switch spikeTypeFlag
    case 's'
        times = s;
    case 'seh'
        
        times = seh;
    case 'd'
        times = d;
    case 'deh'
        times = deh;
    case 't'
        times = t;
    case 'teh'
        times = teh;
    case 'q'
        times = q;
    case 'qeh'
        times = qeh;
end
win = [times-winSize;times+winSize]';
if size(Xin,2)>size(Xin,1)
    Xin = Xin';
end
% 
% fprintf('Scaling each input vector\n')
% for ii = 1:size(Xin,2)
%     Xin(:,ii) = scale(Xin(:,ii));
% end
% 

dum = [];
for ii = 1:size(Xin,2)
    for jj = 1:length(times)
        dum(jj,:) = Xin(times(jj)-winSize:times(jj)+winSize,ii);
    end
    sta(:,ii) = nanmean(dum);
    se(:,ii) = nanstd(dum)./sqrt(length(times));
end
for ii = 1:size(sta,2)
    plot(sta(:,ii));
    ho
%     shadedErrorBar(1:size(sta,1),sta(:,ii),se(:,ii))
end


vline(winSize+1);

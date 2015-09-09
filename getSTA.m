function sta = getSTA(matrix,winSize,varargin)
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



mechTitle = {'FX','FY','M'};
geoTitle = {'R','\Theta'};
if ~ismember(spikeTypeFlag,{'s','seh','d','deh','t','teh','q','qeh'})
    spikeTypeFlag = 's';
    disp('automatically using singlets')
end

s = strfind(matrix(:,1)',1);
seh = strfind(matrix(:,1)',[0 1 0]);seh = seh +1;
d =strfind(matrix(:,1)',[1 1]);
deh = strfind(matrix(:,1)',[0 1 1 0]);deh = deh+1;
t = strfind(matrix(:,1)',[1 1 1]);
teh = strfind(matrix(:,1)',[0 1 1 1 0]);teh = teh+1;
q = strfind(matrix(:,1)',[1 1 1 1]);
qeh = strfind(matrix(:,1)',[0 1 1 1 1 0]);qeh = qeh+1;

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

for ii = 1:length(times)
    trigFX(ii,:) = matrix(times(ii)-winSize:times(ii)+winSize,2);
    trigFY(ii,:) = matrix(times(ii)-winSize:times(ii)+winSize,3);
    trigM(ii,:) = matrix(times(ii)-winSize:times(ii)+winSize,4);
end
se_trigFX = nanstd(trigFX)/sqrt(nansum(matrix(:,1)));
se_trigFY = nanstd(trigFY)/sqrt(nansum(matrix(:,1)));
se_trigM = nanstd(trigM)/sqrt(nansum(matrix(:,1)));

trigFX = nanmean(trigFX);
trigFY = nanmean(trigFY);
trigM = nanmean(trigM);

sta.fx = trigFX;
sta.fx_se = se_trigFX;

sta.fy = trigFY;
sta.fy_se = se_trigFY;

sta.m = trigM;
sta.m_se = se_trigM;


%  
% for ii = 1:size(matrix,2)-1
%     subplot(1,size(matrix,2)-1,ii)
%     
%     
%     plot(-timeLag:-1,staUD(:,ii),'-o')
%     xlabel('ms','FontSize',16)
%     ylabel('newtons','FontSize',16)
%     ln3
%     ho
%     errorbar(-timeLag:-1,staUD(:,ii),staSTD(:,ii));
%     if strcmp(PlotFlag,'mech')
%         title(mechTitle{ii})
%     elseif strcmp(PlotFlag,'geo')
%         title(geoTitle{ii})
%     end
%     
%     
% end
% 

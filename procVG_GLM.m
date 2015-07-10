function out = procVG_GLM(in)
%% parameters
numK = 5;
%%
if size(in.X,2)>size(in.X,1)
    in.X = in.X';
end

if isrow(in.spikes);in.spikes = in.spikes'; end
if isrow(in.C);in.C = in.C'; end


%set first and last times to 0
in.X(1,:) = 0;
in.X(end,:)=0;

% find indices of contact periods
in.C(1) = 0;
in.C(end) = 0;
starts = find(diff(in.C)==1)+1;
stops = find(diff(in.C)==-1);

%set noncontact vars to zeros
in.X(~in.C,:)=0;
in.X_noHist(~in.C,:) = 0;
% interpolate over remaining nans
in.X = naninterp(in.X);
in.X_noHist = naninterp(in.X_noHist);

% scale X
in.X = zscore(in.X);
in.X_noHist = zscore(in.X_noHist);

% init crossvalidate
train.k = crossvalind('Kfold',length(in.C),numK);

% create testing set
test.bool =false(size(in.C));

%this is dumb. Maybe refactor this later
for ii = 1:length(starts)
    
    if starts(ii)<=in.filtSize
        if stops(ii)+in.filtSize>length(test.bool)
            test.bool(1:end) = 1;
            continue
        else
            test.bool(1:stops(ii)+in.filtSize) = 1;
            continue
        end
    end
    if stops(ii)+in.filtSize>length(test.bool)
        test.bool(starts(ii)-in.filtSize:end) = 1;
        continue
    end
    
    test.bool(starts(ii)-in.filtSize:stops(ii)+in.filtSize) = 1;
    
end


test.X = in.X(test.bool,:);
test.X_noHist = in.X_noHist(test.bool,:);
test.k = train.k(test.bool);

%init output
out.Y = zeros(sum(test.bool),1);
out.predictiveModeY = zeros(sum(test.bool),1);

% run GLM
for ii = 1:numK
    [w,dev,stats] = glmfit(in.X(train.k~=ii,:),in.spikes(train.k~=ii),'binomial');
    %wHist = glmfit(in.X(train.k~=ii,:),in.spikes(train.k~=ii),'binomial');
    
    % currently only fitting once . may need to split the history and
    % stimulus fitting as I was doing before
    weights = buildGLM.combineWeights(in.dm,w(2:end));
    if any(strcmp(fieldnames(weights),'X'))
        stimWeights = weights.X.data;
    else
        stimWeights = [];
    end
    if any(strcmp(fieldnames(weights),'derivative'))
        derivWeights = weights.derivative.data;
    else
        derivWeights = [];
    end
    allStimWeight{ii} = stimWeights;
    allDerivWeight{ii} = derivWeights;
    
    allDev{ii} = dev;
    allStats{ii} = stats;
    
    
    out.predictiveModeY(test.k == ii) = glmval(w,test.X(test.k==ii,:),'logit');
    out.Y(test.k == ii) = glmval(w([1 in.histSize+2:end]),test.X(test.k==ii,(in.histSize+1:end)),'identity');
end


% sim GLM
% calculating spike history term on the whole trace
wHist= glmfit(in.X,in.spikes,'binomial');
hist = buildGLM.combineWeights(in.dm,wHist(2:end));
hist = hist.hist.data;

out.raster = simGLM4(out.Y,hist,500);

out.P = mean(out.raster,2);
[x,y,thresh,AUC] = perfcurve(in.spikes(test.bool),out.P,'1');

%% output handling
out.AUC = AUC;
out.ROC_x =x;
out.ROC_y = y;
out.thresh = thresh;
out.stimWeights = allStimWeight;
out.derivWeight = allDerivWeight;
out.dev = allDev;
out.stats = allStats;
out.tested_spikes=  in.spikes(test.bool);
out.spikeHistory = hist;
out.test_set = test.bool;











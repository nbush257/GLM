function out = procVG_GLM(in,bases)
%% parameters
numK = 5;
numHistDims = in.histSize;
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
out.Y_noHist = zeros(sum(test.bool),1);
out.predictiveModeY = zeros(sum(test.bool),1);

% run GLM
for ii = 1:numK
    wNoHist = glmfit(in.X_noHist(train.k~=ii,:),in.spikes(train.k~=ii),'binomial');
    wHist = glmfit(in.X(train.k~=ii,:),in.spikes(train.k~=ii),'binomial');
    
    weightsNoHist = buildGLM.combineWeights(in.dm_noHist,wNoHist(2:end));
    weightsHist = buildGLM.combineWeights(in.dm,wHist(2:end));
    
    if any(strcmp(fieldnames(weightsHist),'X'))
        stimWeights = weightsHist.X.data;
        stimWeightsNoHist = weightsNoHist.X.data;
    else
        stimWeights = [];
        stimWeightsNoHist = [];
    end
    if any(strcmp(fieldnames(weightsHist),'derivative'))
        derivWeights = weightsHist.derivative.data;
        derivWeightsNoHist = weightsNoHist.derivative.data;
    else
        derivWeights = [];
        derivWeightsNoHist = [];
    end
    
    
    allStimWeight{ii} = stimWeights;
    allStimWeightNoHist{ii} = stimWeightsNoHist;
    
    allDerivWeight{ii} = derivWeights;
    allDerivWeightNoHist{ii} = derivWeightsNoHist;
    
    
    out.predictiveModeY(test.k == ii) = glmval(wHist,test.X(test.k==ii,:),'logit');
    out.Y_noHist(test.k ==ii) = glmval(wNoHist,test.X_noHist(test.k==ii,:),'logit');
    out.Y(test.k == ii) = glmval(wHist([1 numHistDims+2:end]),test.X_noHist(test.k==ii,:),'identity');
end


% sim GLM
% calculating spike history term on the whole trace


hist = buildGLM.combineWeights(in.dm,wHist(2:end));
hist = hist.hist.data;

out.raster = simGLM4(out.Y,hist,500);

out.P = mean(out.raster,2);
try
    [x,y,thresh,AUC] = perfcurve(in.spikes(test.bool),out.P,'1');
catch
    x = [];
    y = [];
    thressh = []
    AUC = [];
end

%% output handling
out.AUC = AUC;
out.ROC_x =x;
out.ROC_y = y;
out.thresh = thresh;
out.stimWeights = allStimWeight;
out.stimWeightsNoHist = allStimWeightNoHist;
out.derivWeight = allDerivWeight;
out.derivWeightNoHist = allDerivWeightNoHist;
out.wHist = wHist;
out.wNoHist = wNoHist;
out.dm = in.dm;
out.dm_nh = in.dm_noHist;
out.tested_spikes=  in.spikes(test.bool);
out.spikeHistory = hist;
out.test_set = test.bool;
out.bases = bases;











function [ws,w,X,stats,dev]= runGLM(data,C,spikes, varargin)
%% function [y,w,ws]= runGLM(data,C,spikes,[deriv],[hist])
deriv = 0;
hist = 0;
winSize = 10;
bSize = 4;
offset = 1;
covar = 1;
kfold = 1;
distr = 'binomial';
%% handle inputs
if ~isempty(varargin),
    for ii = 1:2:(length(varargin))
        switch varargin{ii},
            case 'deriv', deriv = varargin{ii+1};
            case 'hist', hist = varargin{ii+1};
            case 'bSize', bSize = varargin{ii+1};
            case 'winSize', winSize = varargin{ii+1};
            case 'offset', offset = varargin{ii+1};
            case 'covar',covar = varargin{ii+1};
            case 'kfold',kfold = varargin{ii+1};
            case 'distr',distr = varargin{ii+1};
            otherwise,
                error('Not a valid input parameter');
        end
    end
end

%%
if size(data,1)<size(data,2)
    data = data';
end

for jj = 1:size(data,2)
    derivative(:,jj) = cdiff(data(:,jj));
end

expt = buildGLM.initExperiment('ms',1,[],[]);
sptrain = find(spikes);
if covar
    expt = buildGLM.registerContinuous(expt,'data',[],size(data,2));
end
if deriv
    expt = buildGLM.registerContinuous(expt,'derivative',[],size(data,2));
end
expt = buildGLM.registerSpikeTrain(expt,'sptrain',[]);

trial = buildGLM.newTrial(expt,length(spikes));
trial.sptrain = sptrain;
if covar
    trial.data = data;
end
if deriv
    trial.derivative = derivative;
end

expt = buildGLM.addTrial(expt,trial,1);
dspec = buildGLM.initDesignSpec(expt);

bsStim = basisFactory.makeSmoothTemporalBasis('raised cosine',winSize, bSize, expt.binfun);

if hist
    dspec = buildGLM.addCovariateSpiketrain(dspec,'hist','sptrain','History Filter');
end

if covar
    dspec = buildGLM.addCovariateRaw(dspec,'data','stimulus',bsStim);
end

if deriv
    dspec = buildGLM.addCovariateRaw(dspec,'derivative','stimulus',bsStim);
end

dm = buildGLM.compileSparseDesignMatrix(dspec,1);

y = buildGLM.getBinnedSpikeTrain(expt,'sptrain',dm.trialIndices);
X = full(dm.X);
y = full(y);
idx = crossvalind('Kfold',length(spikes),kfold);
if kfold>1
    xROC = [];
    yROC = [];
    for kk = 1:kfold
        trainX = X(idx~=kk,:);
        testX = X(idx==kk,:);
        trainY = y(idx~=kk);
        testY = y(idx==kk);
        if strcmp(distr,'binomial')
            [w{kk}, dev{kk}, stats{kk}] = glmfit(trainX, trainY, 'binomial', 'link', 'logit');
            w_noC = w{kk}(2:end);
            ws{kk} = buildGLM.combineWeights(dm,w_noC);
            yhat{kk} = glmval(w{kk},testX,'logit');
            [xROC,yROC,~,AUC(kk)] = perfcurve(testY,yhat{kk},'1');
            plot(xROC,yROC);ho;
        elseif strcmp(distr,'poisson')
            [w{kk}, dev{kk}, stats{kk}] = glmfit(trainX, trainY, 'binomial', 'link', 'logit');
            yhat{kk} = glmval(w{kk},testX,'log');
            corrcoef(testY,yhat{kk});
            r(kk) = corrcoef(1,2);
            
        else
            error('Bad distribution input')
        end
    end
    if strcmp(distr,'binomial')
        [~,keeper] = max(AUC);
        yhat = glmval(w{keeper},dm.X,'logit');
        w = w{keeper};
        ws = buildGLM.combineWeights(dm,w(2:end));
    else
        [~,keeper] = max(r);
        yhat = glmval(w{keeper},dm.X,'log');
        w = w{keeper};
        ws = buildGLM.combineWeights(dm,w(2:end));
    end
    
else
    if strcmp(distr,'poisson')
        [w, dev, stats] = glmfit(X, y, 'poisson', 'link', 'log');
        w_noC = w(2:end);
        ws = buildGLM.combineWeights(dm,w_noC);
%         yhat = glmval(w,dm.X,'logit');
    elseif strcmp(distr,'binomial')
        
        [w, dev, stats] = glmfit(X, y, 'binomial', 'link', 'logit');
        w_noC = w(2:end);
        ws = buildGLM.combineWeights(dm,w_noC);
%         yhat = glmval(w,dm.X,'logit');
    else
        error('Bad distribution input')
    end
end


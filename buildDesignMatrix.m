function [Xout,dm] = buildDesignMatrix(X,spikes,varargin)
%warning('Deirvatives are being heavily smoothed, if you don''t want that you need to edit the deriavtive section')
deriv = 1;
hist = 0;
lastPeak = 50;
bSize = 4;
offset = 1;
covar = 1;
nlOffset = 10;
spikePeak1= 1;
spikePeak2 = 3;

bsStim = basisFactory.makeNonlinearRaisedCos(bSize,1,[offset lastPeak],nlOffset);
bsSpike = basisFactory.makeNonlinearRaisedCos(2,1,[spikePeak1 spikePeak2],nlOffset);
%% handle inputs
if ~isempty(varargin),
    for ii = 1:2:(length(varargin))
        switch varargin{ii},
            case 'deriv', deriv = varargin{ii+1};
            case 'hist', hist = varargin{ii+1};
            case 'bSize', bSize = varargin{ii+1};
            case 'lastPeak', lastPeak = varargin{ii+1};
            case 'offset', offset = varargin{ii+1};
            case 'covar',covar = varargin{ii+1};
            case 'nlOffset',nlOffset = varargin{ii+1};
            case 'bsStim',bsStim = varargin{ii+1};
            case 'bsSpike',bsSpike = varargin{ii+1};
            otherwise,
                error('Not a valid input parameter');
        end
    end
end

%%
if size(X,1)<size(X,2)
    X = X';
end

for jj = 1:size(X,2)
    derivative(:,jj) = cdiff(X(:,jj));
    if any(isnan(derivative(:,jj)))
        warning('NaN in mechanical trace. Derivatives are useless. Turning off derivatives')
        deriv = 0;
    end
        
    derivative(:,jj) = bwfilt(derivative(:,jj),1000,0,10);
end

expt = buildGLM.initExperiment('ms',1,[],[]);
sptrain = find(spikes);
if covar
    expt = buildGLM.registerContinuous(expt,'X',[],size(X,2));
end
if deriv
    expt = buildGLM.registerContinuous(expt,'derivative',[],size(X,2));
end
expt = buildGLM.registerSpikeTrain(expt,'sptrain',[]);

trial = buildGLM.newTrial(expt,length(spikes));
trial.sptrain = sptrain;
if covar
    trial.X = X;
end
if deriv
    trial.derivative = derivative;
end

expt = buildGLM.addTrial(expt,trial,1);
dspec = buildGLM.initDesignSpec(expt);

if hist
    
    dspec = buildGLM.addCovariateSpiketrain(dspec,'hist','sptrain','History Filter',bsSpike);
end

if covar
    dspec = buildGLM.addCovariateRaw(dspec,'X','stimulus',bsStim);
end

if deriv
    dspec = buildGLM.addCovariateRaw(dspec,'derivative','stimulus',bsStim);
end

dm = buildGLM.compileSparseDesignMatrix(dspec,1);

Xout = full(dm.X);

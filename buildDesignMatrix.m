function [Xout,dm] = buildDesignMatrix(X,spikes,varargin)

deriv = 1;
hist = 0;
lastPeak = 50;
bSize = 4;
offset = 1;
covar = 1;
nlOffset = 10;

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

bsStim = basisFactory.makeNonlinearRaisedCos(bSize,1,[offset lastPeak],nlOffset);

if hist
    dspec = buildGLM.addCovariateSpiketrain(dspec,'hist','sptrain','History Filter');
end

if covar
    dspec = buildGLM.addCovariateRaw(dspec,'X','stimulus',bsStim);
end

if deriv
    dspec = buildGLM.addCovariateRaw(dspec,'derivative','stimulus',bsStim);
end

dm = buildGLM.compileSparseDesignMatrix(dspec,1);

Xout = full(dm.X);

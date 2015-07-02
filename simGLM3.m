function [history,lambda] = simGLM3(x,historyFilter,numReplicates)
% input is an 'indentity' link function'
historyFilter = flipud(historyFilter);
plotTGL =1;
xR = repmat(x,1,numReplicates);
spikeOut = zeros(size(xR));
lambda = zeros(size(x));
history = zeros(size(x));
for ii = 2:length(x)
    if ii<=length(historyFilter)
        s = history(1:ii-1);
        lambda(ii) = exp(x(ii) + historyFilter(end-ii+2:end)'*s)./(1+exp(x(ii) + historyFilter(end-ii+2:end)'*s));
    else
        s = history(ii-length(historyFilter):ii-1);
        lambda(ii) = exp(x(ii)+historyFilter'*s)./(1+exp(x(ii)+historyFilter'*s));
    end
    history(ii) = binornd(1,lambda(ii));
end


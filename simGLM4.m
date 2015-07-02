function [history,lambda] = simGLM4(x,historyFilter,numReplicates)
% input is an 'indentity' link function'
historyFilter = flipud(historyFilter);
historyFilter = repmat(historyFilter,1,numReplicates);
plotTGL =0;
xR = repmat(x,1,numReplicates);
lambda = zeros(length(x),numReplicates);
history = zeros(length(x),numReplicates);
for ii = 2:length(x)
    if ii<=size(historyFilter,1)
        s = history(1:ii-1,:);
        lambda(ii,:) = exp(xR(ii,:) + bsxfun(@dot,historyFilter(end-ii+2:end,:),s))./(1+exp(x(ii) + bsxfun(@dot,historyFilter(end-ii+2:end,:),s)));
    else
        s = history(ii-size(historyFilter,1):ii-1,:);
        lambda(ii,:) = exp(xR(ii,:)+bsxfun(@dot,historyFilter,s))./(1+exp(xR(ii,:)+bsxfun(@dot,historyFilter,s)));
    end
    history(ii,:) = binornd(1,lambda(ii),1,numReplicates);
end

if plotTGL
    figure
    for ii = 1:numReplicates
        plot(find(spikeOut(:,ii)==1),ones(size(find(spikeOut(:,ii)==1)))*ii,'k.');ho;
    end
end

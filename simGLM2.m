function spikeOut = simGLM2(x,hist,numReplicates)

plotTGL =1;
xR = repmat(x,1,numReplicates);
spikeOut = zeros(size(xR));
for ii =1:size(xR,1)
    if xR(ii,1)==0
        continue
        if xR(ii,1)>1
            xR(ii,:) = 1;
        end
    end
    %     spike = poissrnd(xR(ii,:));
    spike = binornd(1,xR(ii,:));
    if ii+length(hist)>size(xR,1)
        xR(ii+1:end,find(spike==1)) = xR(ii+1:end,find(spike==1)).*repmat(hist(1:size(xR,1)-ii),1,length(find(spike==1)));
    else
        xR(ii+1:ii+length(hist),find(spike==1)) = xR(ii+1:ii+length(hist),find(spike==1)).*repmat(hist,1,length(find(spike==1)));
    end
    spikeOut(ii,:) = spike;
end

if plotTGL
    figure
    for ii = 1:numReplicates
        plot(find(spikeOut(:,ii)==1),ones(size(find(spikeOut(:,ii)==1)))*ii,'k.');ho;
    end
end

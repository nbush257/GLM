function spikeOut = simGLM(x,hist,numReplicates)
plotTGL = 1;
spikeOut = zeros(length(x),numReplicates);

for nn = 1:numReplicates
    yhatIn = x;
    for ii = 1:length(yhatIn)
        spike = poissrnd(yhatIn(ii));
        if spike
            spikeOut(ii,nn)=1;
            if ii+length(hist)>length(yhatIn)
                temphist = hist(1:length(yhatIn)-ii-1);
            else
                temphist =hist;
            end
            yhatIn(ii+1:ii+length(temphist)) = yhatIn(ii+1:ii+length(temphist)).*temphist;
        end
    end
end

if plotTGL
    figure
    for ii = 1:numReplicates
        plot(find(spikeOut(:,ii)),ones(size(find(spikeOut(:,ii))))*ii,'k.');ho;
    end
end
axy(0,numReplicates+3)

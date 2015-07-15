function xOut = gaussianSpikeSmooth(x,sigma)
m = memory;
if length(x)*length(find(x))*8>m.MemAvailableAllArrays*.1
    warning('Splitting to avoid memory overload this might take longer')
    splits = linspace(1,length(x),5);
    xOut = [];
    for ii = 1:length(splits)-1
        
        clear tempXOut
        if ii == length(splits)-1
            chunk = splits(ii):splits(ii+1)-1;
        else
            chunk = splits(ii):splits(ii+1);
        end
        tempTDiff = bsxfun(@minus,chunk,find(x));
        for jj = 1:length(sigma)
            tempXOut(:,jj) = sum(exp(-tempTDiff.^2/(2*sigma(jj)^2)))/sqrt(2*pi*sigma(jj)^2);
        end
        clear tempTDiff
        xOut = [xOut; tempXOut];
    end
else
    tDiff= bsxfun(@minus,1:length(x),find(x));
    xOut = zeros(length(x),length(sigma));
    for ii = 1:length(sigma)
        xOut(:,ii)= sum(exp(-tDiff.^2/(2*sigma(ii)^2)))/sqrt(2*pi*sigma(ii)^2);
    end
end
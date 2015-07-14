function xOut = gaussianSpikeSmooth(x,sigma)
tDiff= bsxfun(@minus,1:length(x),find(x));
xOut= sum(exp(-tDiff.^2/(2*sigma^2)))/sqrt(2*pi*sigma^2);

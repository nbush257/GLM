% get bayesian contours
function bayesPlots(X,y)
% X is a matrix of inputs
% y is a vector of spikes
% P(A|B) = P(B|A)*P(A)/P(B)
sTimes = find(y);
maxLatency = 100;

dX = diff(X);
triggered  = nan(length(sTimes),size(X,2),10);
for ii = 1:length(sTimes)
    for latency = 1:maxLatency
        triggered(ii,:,latency) = X(sTimes(ii)-latency,:);
    end
end

for ii = 1:maxLatency
    [n,c] = hist3(triggered(:,2:3,ii),[100 100]);
    
    image(n)
    colorbar
    pause(.2)
end
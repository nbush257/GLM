% get bayesian contours
function bayesPlots(X,y)
% X is a matrix of inputs
% y is a vector of spikes
% P(A|B) = P(B|A)*P(A)/P(B)
latencies = [5];
toRM = abs(X)>repmat(std(X)*10,length(y),1);
X(toRM) = 0;
toRM = all(X==0,2);
X(toRM,:) = [];
y(toRM) = [];

sTimes = find(y);
dX = diff(X);
triggered  = nan(length(sTimes),size(X,2),10);
not_triggered = {};
for latency= latencies
    not_triggered{latency} = [zeros(latency,size(X,2)); X(1:end-latency,:)];
end



for ii = 1:length(sTimes)
    for latency = latencies
        triggered(ii,:,latency) = X(sTimes(ii)-latency,:);
        
    end
end

for ii =latencies
    figure(1)
    subplot(221)
    
    ax1 = gca;
    
    [nt,ct] = hist3(triggered(:,1:2,ii),[200 200]);
    nt = nt./length(sTimes)
    
    imagesc(ct{1},ct{2},nt)
    title('Spike Triggered')
     xlabel('\theta')
    ylabel('\phi')
    
    colorbar
    subplot(222)
    
    
    ax2 = gca
    
   
    [nn,cn] = hist3(not_triggered{ii}(:,1:2),[200 200]);
    nn = nn./(length(y));
    imagesc(cn{1},cn{2},nn)
    title('Prior (Stimulus Distribution)')
     xlabel('\theta')
    ylabel('\phi')
    colorbar
    n = nt./nn.*(length(sTimes)/length(y));
    
    subplot(2,2,3)
    ax3 = gca;
    
    imagesc(ct{1},ct{2},n)
    title('Conditional Probability')
     xlabel('\theta')
    ylabel('\phi')
    
%     linkaxes([ax1,ax2,ax3],'xy')
    
    colormap jet
    colorbar
    pause(.2)
    
end
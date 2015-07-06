function kout = kfoldWhisk(C,K,filtSize)
if isrow(C);C = C';end
% make cross validations over contact periods.
start = find(diff(C)==1)+1;
stop = find(diff(C)==-1);
kout = NaN(size(C));
if length(start)<2*K
    warning('number of contacts is less than twice the number of groups. Switching to MATLABs crossvalind')
    kout = crossvalind('Kfold',length(C),K);
else
    for ii = 1:length(start)
        k = randi(K);
        kout(start(ii)-filtSize-5:stop(ii)+filtSize+5) = k;
    end
    
    kout(isnan(kout)') = randi(K, 1,length(kout(isnan(kout))));
end

d = dir('*merged.mat')
for dd = 1:length(d)
    clearvars -except d dd
    load(d(dd).name)
    hist(geo_85(1,:),600)
    [x,~] = ginput(2);
    
    
    C(isnan(C)) = 0;
    C = logical(C);
    blipsUp = strfind(C,[0 1 0])+1;
    blipsDown = strfind(C,[1 0 1])+1;
    C(blipsUp) = 0;
    C(blipsDown) = 1;
    spikevec(spikevec>1)=1;
    
    
    
    TH = geo_85(2,:);
    
    start = find(diff(C)==1)+1;
    stop = find(diff(C)==-1);
    
    
    THnew = nan(size(TH));
    for ii = 1:length(start)
        non_nan = find(~isnan(TH(start(ii):stop(ii))));
        if isempty(non_nan)
            continue
        end
        THnew(start(ii):stop(ii)) = TH(start(ii):stop(ii))-TH(start(ii)+non_nan(1)-1);
    end
    geo_85(2,:) = THnew;
    
    if isempty(x)
        save([d(dd).name(1:end-4) '_cleaned.mat'])
        
        continue
    else
        
        prox = logical(zeros(length(C),1));
        med = logical(zeros(length(C),1));
        dis = logical(zeros(length(C),1));
        
        prox(geo_85(1,:)<x(1))=1;M
        dis(geo_85(1,:)>x(2))=1;
        med(geo_85(1,:)>x(1) & geo_85(1,:)<x(2))=1;
        save([d(dd).name(1:end-4) '_cleaned.mat'])
        
    end
end


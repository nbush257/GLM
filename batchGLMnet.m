d = dir('*toGLM.mat');
for ddd = 1:length(d)
    %% load and clean data
    clearvars -except d ddd
    preRate = 0;
    noProx = 0;
    saveTGL =1;
    d = dir('*toGLM.mat');
    filtSize = 25;
    basisSize = 4;
    rateBin = 15;
    histSize = 15;
    d = dir('*toGLM.mat');
    ca
    fname = [d(ddd).name(1:end-10) '_simGLM'];
    load(d(ddd).name,'mech*','geo*','C','spike*','med','prox','dis')
    if size(mech_85,2)>size(mech_85,1)
        mech_85 = mech_85';
        geo_85 = geo_85';
    end
    
    if ~exist('prox','var')
        prox = zeros(size(C));
    end
    
    if isrow(spikevec);spikevec = spikevec'; end;
    if isrow(C);C = C'; end;
    
    
    if exist('med','var')
        if isrow(med);med = med';end
        if isrow(dis);dis = dis';end
        if isrow(prox);prox = prox';end
    end
    if noProx
        C(prox==1)=0;
    end
    mech_85(1,:) = 0;mech_85(end,:) = 0; mech_85 = naninterp(mech_85);
    geo_85(1,:) = 0;geo_85(end,:) = 0; geo_85 = naninterp(geo_85);
    C(1) = 0;C(end)=0;
    starts = find(diff(C)==1)+1;
    stops = find(diff(C)==-1);
    newMech = [];newGeo = [];newC = [];newSpikes = [];newDis = [];newMed = [];newProx = [];
    for ii = 1:length(starts)
        newMech = [newMech;zeros(filtSize+5,3);mech_85(starts(ii):stops(ii),:);zeros(filtSize+5,3)];
        newGeo = [newGeo;zeros(filtSize+5,2);geo_85(starts(ii):stops(ii),:);zeros(filtSize+5,2)];
        newC = [newC;zeros(filtSize+5,1);C(starts(ii):stops(ii));zeros(filtSize+5,1)];
        newSpikes = [newSpikes;zeros(filtSize+5,1);spikevec(starts(ii):stops(ii));zeros(filtSize+5,1)];
        if exist('med','var')
            newDis = [newDis;ones(filtSize+5,1)*dis(starts(ii)); dis(starts(ii):stops(ii));ones(filtSize+5,1)*dis(stops(ii))];
            newMed = [newMed;ones(filtSize+5,1)*med(starts(ii)); med(starts(ii):stops(ii));ones(filtSize+5,1)*med(stops(ii))];
            newProx = [newProx;ones(filtSize+5,1)*prox(starts(ii)); prox(starts(ii):stops(ii));ones(filtSize+5,1)*prox(stops(ii))];
            newDis = logical(newDis);newMed = logical(newMed);newProx = logical(newMed);
        end
    end
    
    
    if sum(newSpikes)<20
        warning('Not enough spikes found in the region of interest. Aborting.')
        continue
    end
    %% calculate rate pre GLM
    if preRate
        rate = tsmovavg(newSpikes','s',rateBin)*rateBin;rate = rate';
        rateMech = tsmovavg(newMech','s',rateBin);rateMech = rateMech';
        rateGeo = tsmovavg(newGeo','s',rateBin);rateGeo = rateGeo';
        rate(isnan(rate)) = 0;rate = rate*1000;
        rateMech(isnan(rateMech))=0;
        rateGeo(isnan(rateGeo))=0;
        
        %in case we want to actually downsample our binsize
        rate = [];rateMech=[];rateGeo = [];
        ii=0;startBin = 1;stopBin=startBin+rateBin-1;
        while stopBin <=length(newSpikes)
            ii = ii+1;
            rate(ii) = sum(newSpikes(startBin:stopBin));
            rateMech(ii,:) = mean(newMech(startBin:stopBin,:));
            rateGeo(ii,:) = mean(newGeo(startBin:stopBin,:));
            startBin = ii*rateBin+1;stopBin = startBin + rateBin-1;
        end
        if isrow(rate);rate = rate';end
        rate = rate*1000;
        
        
        [XMr,dmMr] = buildDesignMatrix(rateMech,rate,'winSize',filtSize,'bSize',basisSize);
        fitMr = cvglmnet(XMr,rate,'poisson',[],[],5,[]);
        
        [XGr,dmGr] = buildDesignMatrix(rateGeo,rate,'winSize',filtSize,'bSize',basisSize);
        fitGr = cvglmnet(XGr,rate,'poisson',[],[],5,[]);
        YMr = cvglmnetPredict(fitMr,XMr,[],'response');
        YGr = cvglmnetPredict(fitGr,XGr,[],'response');
        
        wMr = cvglmnetCoef(fitMr);
        wGr = cvglmnetCoef(fitGr);
        weightMr = buildGLM.combineWeights(dmMr,wMr(2:end));
        weightGr= buildGLM.combineWeights(dmGr,wGr(2:end));
        rGr = corr(YGr,rate);
        rMr = corr(YMr,rate);
        
    end
    
    %% fit GLM
    
    k = crossvalind('Kfold',length(newSpikes),5);
    [XM,dmM] = buildDesignMatrix(newMech,newSpikes,'winSize',filtSize,'bSize',basisSize,'hist',0);
    [XG,dmG] = buildDesignMatrix(newGeo,newSpikes,'winSize',filtSize,'bSize',basisSize,'hist',0);
    [XMh,dmMh] = buildDesignMatrix(newMech,newSpikes,'winSize',filtSize,'bSize',basisSize,'hist',1);
    [XGh,dmGh] = buildDesignMatrix(newGeo,newSpikes,'winSize',filtSize,'bSize',basisSize,'hist',1);

    for ii = 1:max(k)+1
        ii
        if ii>max(k)
            
            [~,keepG] = max(rG);wG = allWG{keepG};wGh = allWGH{keepG};
            [~,keepM] = max(rM);wM = allWM{keepM};wMh = allWMH{keepM};
            k = ones(size(newSpikes))*ii;
        else
            wMh = glmfit(XMh(k~=ii,:),newSpikes(k~=ii),'binomial');
            wM = glmfit(XM(k~=ii,:),newSpikes(k~=ii),'binomial');
            wG =glmfit(XG(k~=ii,:),newSpikes(k~=ii),'binomial');
            wGh = glmfit(XGh(k~=ii,:),newSpikes(k~=ii),'binomial');
            
        end
        histM = buildGLM.combineWeights(dmMh,wMh(2:end));histM = histM.hist.data;
        histG = buildGLM.combineWeights(dmGh,wGh(2:end));histG = histG.hist.data;
        
        weightM = buildGLM.combineWeights(dmM,wM(2:end));
        weightG = buildGLM.combineWeights(dmG,wG(2:end));
        
        YM = glmval(wM,XM(k==ii,:),'identity');
        YG = glmval(wG,XG(k==ii,:),'identity');
        %% sim trials
        
        mechOut = simGLM4(YM,histM(1:histSize),500);
        geoOut = simGLM4(YG,histG(1:histSize),500);
        
        mechRate = tsmovavg(mechOut','s',rateBin);mechRate = nanmean(mechRate);mechRate = mechRate';mechRate(isnan(mechRate))=0;mechRate = mechRate*1000;
        geoRate = tsmovavg(geoOut','s',rateBin);geoRate = nanmean(geoRate);geoRate = geoRate';geoRate(isnan(geoRate))=0;geoRate = geoRate*1000;
        rate = tsmovavg(newSpikes','s',rateBin);rate = rate';rate(isnan(rate))=0;rate = rate*1000;
        rate = rate(k==ii);
        
        
        hM = histM(1:histSize);
        hG = histG(1:histSize);
        
        %% get correlations
        rG(ii) = corr(geoRate,rate);
        rM(ii) = corr(mechRate,rate);
        if exist('med','var')
            
            rM_dis = corr(mechRate(newDis),rate(newDis));
            rM_med = corr(mechRate(newMed),rate(newMed));
            
            rG_dis = corr(geoRate(newDis),rate(newDis));
            rG_med = corr(geoRate(newMed),rate(newMed));
            if ~noProx
                rG_prox = corr(geoRate(newProx),rate(newProx));
                rM_prox = corr(mechRate(newProx),rate(newProx));
            end
        end
        
        allWMH{ii} = wMh;
        allWM{ii} = wM;
        allWGH{ii} = wGh;
        allWG{ii} = wG;
    end
    rG = rG(end);
    rM = rM(end);
    weightM = buildGLM.combineWeights(dmM,wM(2:end));
    weightG = buildGLM.combineWeights(dmG,wG(2:end));
    
    histM = buildGLM.combineWeights(dmMh,wMh(2:end));histM = histM.hist.data;
    histG = buildGLM.combineWeights(dmGh,wGh(2:end));histG = histG.hist.data;
    
    %% plot
    f1 = figure;
    subplot(211)
    plot(geoRate);ho;plot(rate);legend({'Predicted','Actual'});
    title(['Geometry: Pearson Correlation Coefficent = ' num2str(rG)])
    xlabel('time (ms)')
    ylabel('Spike Rate')
    subplot(212)
    
    plot(mechRate);ho;plot(rate);legend({'Predicted','Actual'});
    title(['Mechanics: Pearson Correlation Coefficent = ' num2str(rM)])
    xlabel('time (ms)')
    ylabel('Spike Rate')
    
    
    if preRate
        f2 = figure;
        
        plot(YGr);ho;plot(rate);legend({'Predicted','Actual'});
        title(['Geometry: Pearson Correlation Coefficent = ' num2str(rGr)])
        xlabel('time (ms)')
        ylabel('rate (spikes/second)]')
        
        
        plot(YMr);ho;plot(rate);legend({'Predicted','Actual'});
        title(['Mechanics: Pearson Correlation Coefficent = ' num2str(rMr)])
        xlabel('time (ms)')
        ylabel('rate (spikes/second)]')
        if saveTGL
            cd done\
            saveas(f2,[fname '_rates.fig'],'fig')
            cd ..
        end
    end
    
    f3 = figure;
    subplot(231)
    plot(weightG.X.data);title('Geometry Filter');legend({'R','\Theta Push'});
    subplot(232)
    plot(weightG.derivative.data);title('Geometry Time Derivative Filter');legend({'R-dot','\Theta Push-dot'});
    subplot(233)
    plot(hG); title('Geometry Spike Filter')
    
    subplot(234)
    plot(weightM.X.data);title('Mechanics Filter');legend({'FX','FY','M'});
    subplot(235)
    plot(weightM.derivative.data);title('Mechanics Time Derivative Filter');legend({'FX-dot','FY-dot','M-dot'});
    subplot(236)
    plot(hM);title('Mechanics spike history');
    
    if exist('med','var')
        
        f4 = figure;
        subplot(232)
        plot(geoRate(find(newMed)));ho; plot(rate(find(newMed)));title(['Geometry Medial, R = ' num2str(rG_med)])
        subplot(233)
        plot(geoRate(find(newDis)));ho; plot(rate(find(newDis)));title(['Geometry Distal, R = ' num2str(rG_dis)])
        if ~noProx
            
            subplot(231)
            plot(geoRate(find(newProx)));ho; plot(rate(find(newProx)));title(['Geometry Proximal, R = ' num2str(rG_prox)])
            subplot(234)
            plot(mechRate(find(newProx)));ho; plot(rate(find(newProx)));title(['Mechanical Prox, R = ' num2str(rM_prox)])
            
        end
        subplot(235)
        plot(mechRate(find(newMed)));ho; plot(rate(find(newMed)));title(['Mechanics Medial, R = ' num2str(rM_med)])
        subplot(236)
        plot(mechRate(find(newDis)));ho; plot(rate(find(newDis)));title(['Mechanics Distal, R = ' num2str(rM_dis)])
        if saveTGL
            cd done\
            hgsave(f4,[fname '_radialDistances.fig'])
            cd ..
        end
    end
    if saveTGL
        cd done\
        save([fname '.mat']);
        hgsave(f1,[fname '_responses.fig']);
        hgsave(f3,[fname '_filters.fig']);
        cd ..
        ca
        pause(.01)
    else
        pause
    end
    
    
    
end




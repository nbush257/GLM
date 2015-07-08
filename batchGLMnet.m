ca
d = dir('*toGLM.mat');
basisSize = 4;
lastPeak = 10;
nlOffset = 1;
rateBin = 25;
histSize = 5;
bas = basisFactory.makeNonlinearRaisedCos(basisSize,1,[1 lastPeak],nlOffset);
plot(bas.B)
saveTGL =1;% 1 if you want to save the outputs
pause
for ddd = 1:length(d)
    %% load and clean data
    ca
    clearvars -except d ddd cellType lastPeak basisSize nlOffset rateBin histSize saveTGL
    preRate = 0; % 1 if you want to calculate the rate as the predictor before fitting the GLM. Not well tested
    noProx = 1;% 1 if you want to omit proximal deflections.
    d = dir('*toGLM.mat');% get the files in the directory that match the pathspec
    % size of window to bin the rates at
    % number of ms in the past to consider the spike history
    d = dir('*toGLM.mat');
    ca
    fname = [d(ddd).name(1:end-10) '_simGLM'];% output filename -no suffix
    try
        load(d(ddd).name,'mech*','geo*','C','spike*','med','prox','dis')
    catch
        warning('FIle not found')
        pause
        continue
    end
    %     Enable this if you want to do something special for RA cells
    %     if cellType(ddd).type == 1
    %         histSize = 5;
    %     end
    
    % make sure the variables are in the proper orientation ( Time along the
    % first dimension
    if size(mech_85,2)>size(mech_85,1)
        mech_85 = mech_85';
        geo_85 = geo_85';
    end
    
    if ~exist('prox','var')
        prox = zeros(size(C));
    end
    % make sure time and spikes are column vectors
    if isrow(spikevec);spikevec = spikevec'; end;
    if isrow(C);C = C'; end;
    
    
    if exist('med','var')
        if isrow(med);med = med';end
        if isrow(dis);dis = dis';end
        if isrow(prox);prox = prox';end
    end
    
    % ignore all proximal contacts
    if noProx
        C(prox==1)=0;
    end
    % interpolate over nans.
    mech_85(1,:) = 0;mech_85(end,:) = 0; mech_85 = naninterp(mech_85);
    geo_85(1,:) = 0;geo_85(end,:) = 0; geo_85 = naninterp(geo_85);
    
    % find indices of contact periods
    C(1) = 0;C(end)=0; % ensures at least one contact period and that the number of stops = the number of starts
    starts = find(diff(C)==1)+1;
    stops = find(diff(C)==-1);
    % create new variables that consist only of the time of contact + a pad
    % 5 ms longer than the size of the input filter.
    newMech = [];newGeo = [];newC = [];newSpikes = [];newDis = [];newMed = [];newProx = [];
    for ii = 1:length(starts)
        newMech = [newMech;zeros(lastPeak+10,3);mech_85(starts(ii):stops(ii),:);zeros(lastPeak+10,3)];
        newGeo = [newGeo;zeros(lastPeak+10,2);geo_85(starts(ii):stops(ii),:);zeros(lastPeak+10,2)];
        newC = [newC;zeros(lastPeak+10,1);C(starts(ii):stops(ii));zeros(lastPeak+10,1)];
        newSpikes = [newSpikes;zeros(lastPeak+10,1);spikevec(starts(ii):stops(ii));zeros(lastPeak+10,1)];
        if exist('med','var')
            newDis = [newDis;ones(lastPeak+10,1)*dis(starts(ii)); dis(starts(ii):stops(ii));ones(lastPeak+10,1)*dis(stops(ii))];
            newMed = [newMed;ones(lastPeak+10,1)*med(starts(ii)); med(starts(ii):stops(ii));ones(lastPeak+10,1)*med(stops(ii))];
            newProx = [newProx;ones(lastPeak+10,1)*prox(starts(ii)); prox(starts(ii):stops(ii));ones(lastPeak+10,1)*prox(stops(ii))];
            newDis = logical(newDis);newMed = logical(newMed);newProx = logical(newMed);
        end
    end
    % Scale the input vectors
    newMech = zscore(newMech);
    newGeo = zscore(newGeo);
    
    % Short circuit if there are too few spikes
    if sum(newSpikes)<20
        warning('Not enough spikes found in the region of interest. Aborting.')
        continue
    end
    %% calculate rate pre GLM
    % This is not well implemented
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
        
        
        [XMr,dmMr] = buildDesignMatrix(rateMech,rate,'lastPeak',lastPeak,'bSize',basisSize);
        fitMr = cvglmnet(XMr,rate,'poisson',[],[],5,[]);
        
        [XGr,dmGr] = buildDesignMatrix(rateGeo,rate,'lastPeak',lastPeak,'bSize',basisSize);
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
    % prepare the cross validation
    numK=5;
    k = crossvalind('Kfold',length(newSpikes),numK);
    %     k = kfoldWhisk(newC,numK,filtSize);
    
    % Apply the filters and set up the dsign matrices. Uses Pillow's glm
    % code
    [XM,dmM] = buildDesignMatrix(newMech,newSpikes,'lastPeak',lastPeak,'bSize',basisSize,'hist',0);% mechanics
    [XG,dmG] = buildDesignMatrix(newGeo,newSpikes,'lastPeak',lastPeak,'bSize',basisSize,'hist',0);% Geo
    [XMh,dmMh] = buildDesignMatrix(newMech,newSpikes,'lastPeak',lastPeak,'bSize',basisSize,'hist',1);% Mech spike history
    [XGh,dmGh] = buildDesignMatrix(newGeo,newSpikes,'lastPeak',lastPeak,'bSize',basisSize,'hist',1);% Geo spike history
    [XB,dmB] = buildDesignMatrix([newMech newGeo],newSpikes,'lastPeak',lastPeak,'bSize',basisSize,'hist',0);% Both
    [XBh,dmBh] = buildDesignMatrix([newMech newGeo],newSpikes,'lastPeak',lastPeak,'bSize',basisSize,'hist',1);% Both spike history
    % initialize the GLMval output.
    YM = zeros(size(newC));
    YB = zeros(size(newC));
    YG = zeros(size(newC));
    
    for ii = 1:numK % Crossvalidation fitting
        ii
        if ii>numK % should be depricated.
            [~,keepG] = max(rG);wG = allWG{keepG};wGh = allWGH{keepG};
            [~,keepB] = max(rB);wB = allWB{keepB};wBh = allWBH{keepB};
            [~,keepM] = max(rM);wM = allWM{keepM};wMh = allWMH{keepM};
            k = ones(size(newSpikes))*ii;
        else % Fit the glm to the training data
            wMh = glmfit(XMh(k~=ii,:),newSpikes(k~=ii),'binomial');
            wBh = glmfit(XBh(k~=ii,:),newSpikes(k~=ii),'binomial');
            wM = glmfit(XM(k~=ii,:),newSpikes(k~=ii),'binomial');
            wG =glmfit(XG(k~=ii,:),newSpikes(k~=ii),'binomial');
            wB =glmfit(XB(k~=ii,:),newSpikes(k~=ii),'binomial');
            wGh = glmfit(XGh(k~=ii,:),newSpikes(k~=ii),'binomial');
            
        end
        % extract the spike history filter
        histM = buildGLM.combineWeights(dmMh,wMh(2:end));histM = histM.hist.data;
        histG = buildGLM.combineWeights(dmGh,wGh(2:end));histG = histG.hist.data;
        histB = buildGLM.combineWeights(dmBh,wBh(2:end));histB = histB.hist.data;
        % save all the spike history filters from the seperate
        % crossvalidations
        allHistM(:,ii) = histM;
        allHistG(:,ii) = histG;
        allHistB(:,ii) = histB;
        
        % extract all the stimulus weights for each raining set
        weightM(ii)= buildGLM.combineWeights(dmM,wM(2:end));
        weightG(ii)= buildGLM.combineWeights(dmG,wG(2:end));
        weightB(ii)= buildGLM.combineWeights(dmB,wB(2:end));
        
        % Calculate the predicted value for the test set. After looping
        % through all k, this vector will be full.
        YM(k==ii) = glmval(wM,XM(k==ii,:),'identity');
        YG(k==ii) = glmval(wG,XG(k==ii,:),'identity');
        YB(k==ii) = glmval(wB,XB(k==ii,:),'identity');
        
    end
    %% sim trials
    
    % use the mean spike history filter as the history filter in the
    % simulations
    histM = mean(allHistM,2);hM = histM(1:histSize);
    histG = mean(allHistG,2);hG = histG(1:histSize);
    histB = mean(allHistB,2);hB = histB(1:histSize);
    
    % execute the simulations using the predetermined history length
    mechOut = simGLM4(YM,histM(1:histSize),500);
    geoOut = simGLM4(YG,histG(1:histSize),500);
    bothOut = simGLM4(YB,histB(1:histSize),500);
    
    % Calcualte the rate for each sim trial
    mechRate = tsmovavg(mechOut','s',rateBin);mechRate= mechRate*1000;
    % calculate the STD of the rates across trials
    mechSE = nanstd(mechRate);
    %get the mean rate cross all the trials
    mechRate = nanmean(mechRate);mechRate = mechRate';mechRate(isnan(mechRate))=0;
    
    % as above.
    geoRate = tsmovavg(geoOut','s',rateBin);geoRate = geoRate*1000;
    geoSE = nanstd(geoRate);
    geoRate = nanmean(geoRate);geoRate = geoRate';geoRate(isnan(geoRate))=0;
    
    bothRate = tsmovavg(bothOut','s',rateBin);bothRate = bothRate*1000;
    bothSE = nanstd(bothRate);
    bothRate = nanmean(bothRate);bothRate = bothRate';bothRate(isnan(bothRate))=0;
    
    % calculate the spike rate.
    rate = tsmovavg(newSpikes','s',rateBin);rate = rate';rate(isnan(rate))=0;rate = rate*1000;
    
    %% get correlations
    rG = corr(geoRate,rate);
    rM = corr(mechRate,rate);
    rB = corr(bothRate,rate);
    if exist('med','var')
        
        rM_dis = corr(mechRate(newDis),rate(newDis));
        rM_med = corr(mechRate(newMed),rate(newMed));
        
        rG_dis = corr(geoRate(newDis),rate(newDis));
        rG_med = corr(geoRate(newMed),rate(newMed));
        
        rB_dis = corr(bothRate(newDis),rate(newDis));
        rB_med = corr(bothRate(newMed),rate(newMed));
        
        if ~noProx
            rG_prox = corr(geoRate(newProx),rate(newProx));
            rM_prox = corr(mechRate(newProx),rate(newProx));
            rB_prox = corr(bothRate(newProx),rate(newProx));
        end
    end
    
    
    %% plot
    
    % figure of the prediction overlaid with the actual
    f1 = figure;
    subplot(311)
    shadedErrorBar(1:length(rate),geoRate,geoSE);ho; plot(geoRate);ho;plot(rate);legend({'Predicted','Actual'});
    title(['Geometry: Pearson Correlation Coefficent = ' num2str(rG)])
    xlabel('time (ms)')
    ylabel('Spike Rate')
    subplot(312)
    
    shadedErrorBar(1:length(mechRate),mechRate,mechSE);ho;plot(rate);legend({'Predicted','Actual'});
    title(['Mechanics: Pearson Correlation Coefficent = ' num2str(rM)])
    xlabel('time (ms)')
    ylabel('Spike Rate')
    
    subplot(313)
    shadedErrorBar(1:length(bothRate),bothRate,bothSE);ho;plot(rate);legend({'Predicted','Actual'});
    title(['Both: Pearson Correlation Coefficent = ' num2str(rB)])
    xlabel('time (ms)')
    ylabel('Spike Rate')
    
    
    
    if preRate % again, probably depricated. Used if the rates were computed as the inputs to the GLM.
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
    %   Save filters, doesn't work with crossvalidated data.
    %     f3 = figure;
    %     subplot(231)
    %     plot(weightG.X.data);title('Geometry Filter');legend({'R','\Theta Push'});
    %     subplot(232)
    %     plot(weightG.derivative.data);title('Geometry Time Derivative Filter');legend({'R-dot','\Theta Push-dot'});
    %     subplot(233)
    %     plot(hG); title('Geometry Spike Filter')
    %
    %     subplot(234)
    %     plot(weightM.X.data);title('Mechanics Filter');legend({'FX','FY','M'});
    %     subplot(235)
    %     plot(weightM.derivative.data);title('Mechanics Time Derivative Filter');legend({'FX-dot','FY-dot','M-dot'});
    %     subplot(236)
    %     plot(hM);title('Mechanics spike history');
    
    % plot the dependency on radial distance if possible.
    if exist('med','var')
        
        f4 = figure;
        subplot(232)
        shadedErrorBar(1:length(geoRate(find(newMed))),geoRate(find(newMed)),geoSE(find(newMed)));ho; plot(rate(find(newMed)));title(['Geometry Medial, R = ' num2str(rG_med)])
        subplot(233)
        shadedErrorBar(1:length(geoRate(find(newDis))),geoRate(find(newDis)),geoSE(find(newDis)));ho; plot(rate(find(newDis)));title(['Geometry Distal, R = ' num2str(rG_dis)])
        if ~noProx
            
            subplot(231)
            shadedErrorBar(1:length(geoRate(find(newProx))),geoRate(find(newProx)),geoSE(find(newProx)));ho; plot(rate(find(newProx)));title(['Geometry Proximal, R = ' num2str(rG_prox)])
            subplot(234)
            shadedErrorBar(1:length(mechRate(find(newProx))),mechRate(find(newProx)),mechSE(find(newProx)));ho; plot(rate(find(newProx)));title(['Mechanical Prox, R = ' num2str(rM_prox)])
            
        end
        subplot(235)
        shadedErrorBar(1:length(mechRate(find(newMed))),mechRate(find(newMed)),mechSE(find(newMed)));ho; plot(rate(find(newMed)));title(['Mechanics Medial, R = ' num2str(rM_med)])
        subplot(236)
        shadedErrorBar(1:length(mechRate(find(newDis))),mechRate(find(newDis)),mechSE(find(newDis)));ho; plot(rate(find(newDis)));title(['Mechanics Distal, R = ' num2str(rM_dis)])
        if saveTGL
            cd done\
            hgsave(f4,[fname '_radialDistances.fig'])
            cd ..
        end
    end
    
    % save the data if desired.
    if saveTGL
        cd done\
        save([fname '.mat']);
        hgsave(f1,[fname '_responses.fig']);
        %        hgsave(f3,[fname '_filters.fig']);
        cd ..
        pause(.01)
    else
        pause
    end
    
    
    
end





dI = dir('*toGLM.mat');
for ddd =1:length(dI)
    clearvars -except dI ddd
    filtSize = 50;
    saveFlag = 1;
    
    dI = dir('*toGLM.mat');
    ca
    fname = [dI(ddd).name(1:end-10) '_simGLM'];
    load(dI(ddd).name)
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
    C(prox==1)=0;
    mech_85(1,:) = 0;mech_85(end,:) = 0; mech_85 = naninterp(mech_85);
    geo_85(1,:) = 0;geo_85(end,:) = 0; geo_85 = naninterp(geo_85);
    
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
        end
    end
    
    if sum(newSpikes)<25
        warning('Not enough spikes found in the region of interest. Aborting.')
        continue
    end
    %[wsM,~,~] = runGLM(newMech,newC,newSpikes,'hist',1,'deriv',1,'covar',1);
    %hM = wsM.hist.data;
    %%histM = hM;
   % histM = exp(hM);
%    histM = exp(hM)./(1+exp(hM));
    
    [wsM,wM,XM] = runGLM(newMech,newC,newSpikes,'hist',0,'deriv',1,'winSize',filtSize,'bSize',4);
    %YM = glmval(wM([1 12:end]),XM(:,11:end),'log');
    YM = glmval(wM,XM,'logit');
    %histM = histM(1:50);
    %pSpikeM = simGLM2(YM,histM,1000);
    
    
    
   % [wsG,~,~] = runGLM(newGeo,newC,newSpikes,'hist',1,'deriv',1,'covar',1);
   % hG = wsG.hist.data;
   % histG = exp(hG);
   
    [wsG,wG,XG,statsG] = runGLM(newGeo,newC,newSpikes,'hist',0,'deriv',1,'winSize',filtSize,'bSize',4);
    %YG = glmval(wG([1 12:end]),XG(:,11:end),'log');
    YG = glmval(wG,XG,'logit');
    %histG = histG(1:2);
   % histG = exp(hG)./(1+exp(hG));
   % pSpikeG = simGLM2(YG,histG,1000);
    %ca;
    f1 = figure;
    
    
    %tM = tsmovavg(pSpikeM','s',25);
    %tG = tsmovavg(pSpikeG','s',25);
    tM = tsmovavg(YM','s',25);tM(isnan(tM))=0;tM = tM*1000;
    tG = tsmovavg(YG','s',25);tG(isnan(tG))=0;tG = tG*1000;
    
%     
%     mean_tM = nanmean(tM,1)*1000;mean_tM(isnan(mean_tM))=0;
%     mean_tG = nanmean(tG,1)*1000;mean_tG(isnan(mean_tG))=0;
    spikeRate = tsmovavg(newSpikes','s',25)*1000;
    spikeRate(isnan(spikeRate))=0;
    
    % calculate oerrelations
    corrcoef(tG,spikeRate);
    rG = ans(1,2);
    corrcoef(tM,spikeRate);
    rM = ans(1,2);
    if exist('med','var')
        corrcoef(tM(find(newDis)),spikeRate(find(newDis)));
        rM_dis = ans(1,2);
        corrcoef(tM(find(newMed)),spikeRate(find(newMed)));
        rM_med = ans(1,2);
        corrcoef(tG(find(newDis)),spikeRate(find(newDis)));
        rG_dis = ans(1,2);
        corrcoef(tG(find(newMed)),spikeRate(find(newMed)));
        rG_med = ans(1,2);
        corrcoef(tG(find(newProx)),spikeRate(find(newProx)));
        rG_prox = ans(1,2);      
        corrcoef(tM(find(newProx)),spikeRate(find(newProx)));
        rM_prox = ans(1,2); 
    end
    
    
    
    
    
    subplot(211)
    plot(tG);ho;plot(spikeRate);legend({'Predicted','Actual'})
    title(['Geometry: Pearson Correlation Coefficent = ' num2str(rG)])
    xlabel('time (ms)')
    ylabel('rate (spikes/second)[25 ms bins]')
    subplot(212)
    plot(tM);ho;plot(spikeRate);legend({'Predicted','Actual'})
    title(['Mechanics: Pearson Correlation Coefficient = ' num2str(rM)])
    xlabel('time (ms)')
    ylabel('rate (spikes/second) [25ms bins]')
    f2 = figure
    
    subplot(231)
    plot(wsG.data.data);title('Geometry Filter');legend({'R','\Theta Push'});
    subplot(232)
    plot(wsG.derivative.data);title('Geometry Time Derivative Filter');legend({'R-dot','\Theta Push-dot'});
    %subplot(233)
    %plot(histG); title('Geometry Spike Filter')
    
    subplot(234)
    plot(wsM.data.data);title('Mechanics Filter');legend({'FX','FY','M'});
    subplot(235)
    plot(wsM.derivative.data);title('Mechanics Time Derivative Filter');legend({'FX-dot','FY-dot','M-dot'});
    %subplot(236)
    %plot(histM);title('Mechanics spike history');
    f3 = figure;
    
    if exist('med','var')
        
        subplot(231)
        plot(tG(find(newProx)));ho; plot(spikeRate(find(newProx)));title(['Geometry Proximal, R = ' num2str(rG_prox)])
        subplot(232)
        plot(tG(find(newMed)));ho; plot(spikeRate(find(newMed)));title(['Geometry Medial, R = ' num2str(rG_med)])
        subplot(233)
        plot(tG(find(newDis)));ho; plot(spikeRate(find(newDis)));title(['Geometry Distal, R = ' num2str(rG_dis)])
        subplot(234)
        plot(tM(find(newProx)));ho; plot(spikeRate(find(newProx)));title(['Mechanical Prox, R = ' num2str(rM_prox)])
        subplot(235)
        plot(tM(find(newMed)));ho; plot(spikeRate(find(newMed)));title(['Mechanics Medial, R = ' num2str(rM_med)])
        subplot(236)
        plot(tM(find(newDis)));ho; plot(spikeRate(find(newDis)));title(['Mechanics Distal, R = ' num2str(rM_dis)])
        
    end
    if saveFlag
        saveas(f1,[fname '_response.fig'],'fig');
        saveas(f2,[fname '_filters.fig'],'fig')
        saveas(f3,[fname '_CompareRadialDis.fig'],'fig')
        save([fname '.mat'])
        
        ca
        pause(.01)
    end
    
    
end
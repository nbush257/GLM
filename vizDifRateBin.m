rateBin = 50;
d = dir('*Hack.mat')
for dd = 1:length(d)
    load(d(dd).name,'mechOut','geoOut','bothOut','newSpikes')
    mechRate = tsmovavg(mechOut','s',rateBin);mechRate = nanmean(mechRate);mechRate = mechRate';mechRate(isnan(mechRate))=0;mechRate = mechRate*1000;
    geoRate = tsmovavg(geoOut','s',rateBin);geoRate = nanmean(geoRate);geoRate = geoRate';geoRate(isnan(geoRate))=0;geoRate = geoRate*1000;
    bothRate = tsmovavg(bothOut','s',rateBin);bothRate = nanmean(bothRate);bothRate = bothRate';bothRate(isnan(bothRate))=0;bothRate = bothRate*1000;
    rate = tsmovavg(newSpikes','s',rateBin);rate = rate';rate(isnan(rate))=0;rate = rate*1000;
    
    rG = corr(geoRate,rate);
    rM = corr(mechRate,rate);
    rB = corr(bothRate,rate);
    
    
    f1 = figure;
    subplot(311)
    plot(geoRate);ho;plot(rate);legend({'Predicted','Actual'});
    title(['Geometry: Pearson Correlation Coefficent = ' num2str(rG)])
    xlabel('time (ms)')
    ylabel('Spike Rate')
    subplot(312)
    
    plot(mechRate);ho;plot(rate);legend({'Predicted','Actual'});
    title(['Mechanics: Pearson Correlation Coefficent = ' num2str(rM)])
    xlabel('time (ms)')
    ylabel('Spike Rate')
    
    subplot(313)
    plot(bothRate);ho;plot(rate);legend({'Predicted','Actual'});
    title(['Both: Pearson Correlation Coefficent = ' num2str(rB)])
    xlabel('time (ms)')
    ylabel('Spike Rate')
    
    
    
    
    pause(.01)
    c;ca;
end
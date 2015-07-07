rateBin = 25;
d = dir('*simGLM.mat')

type = [0;2;0;2;0;1;2;0;1;1;0;0;1;2;0;2;0;0;0;0;0;1;1;1;0;0;1;1;2;0;0;0;0;0;0;1;1;1;0;0];
strType = {'S','R','M'};
allPrefDir =[NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,90,225,135,270,45,135,NaN,NaN,NaN,135,180,NaN,NaN,NaN,315,180,135,270,225,270,90,NaN,NaN,180,135];
for dd = 1:length(d)
    d(dd).name
    load(d(dd).name,'mechOut','geoOut','bothOut','newSpikes');
    mechRate = tsmovavg(mechOut','s',rateBin);
    mechSE = nanstd(mechRate)./sqrt(size(mechRate,1));
    mechRate = nanmean(mechRate);mechRate = mechRate';mechRate(isnan(mechRate))=0;mechRate = mechRate*1000;
    geoRate = tsmovavg(geoOut','s',rateBin);
    
    geoSE = nanstd(geoRate)./sqrt(size(geoRate,1));
    geoRate = nanmean(geoRate);geoRate = geoRate';geoRate(isnan(geoRate))=0;geoRate = geoRate*1000;
    bothRate = tsmovavg(bothOut','s',rateBin);
    
    bothSE = nanstd(bothRate)./sqrt(size(bothRate,1));
    bothRate = nanmean(bothRate);bothRate = bothRate';bothRate(isnan(bothRate))=0;bothRate = bothRate*1000;
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
    uicontrol('Style','text','String',['Type = ' strType(type(dd)+1)],'Position',[20 0 200 30])
    
    set(gcf,'Position',[ 1          41        1920         963])
    
    pause
    c;ca;
end
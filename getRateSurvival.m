di = dir('*simGLM.mat')

for dIter = 1:length(di)
    clearvars -except di dIter
    
    load(di(dIter).name)
    ca
    rrM = zeros(1,100);
    rrG = zeros(1,100);
    for ii = 1:100
        ntM = tsmovavg(pSpikeM','s',ii);
        ntG = tsmovavg(pSpikeG','s',ii);
        rater = tsmovavg(newSpikes','s',ii);
        ntM = nanmean(ntM,1);
        ntG = nanmean(ntG,1);
        ntM(isnan(ntM))=0;
        ntG(isnan(ntG))=0;
        rater(isnan(rater))=0;
        corrcoef(ntM,rater);
        rrM(ii) = ans(1,2);
        corrcoef(ntG,rater);
        rrG(ii) = ans(1,2);
    end
    %     ca
    %     plot(rrM);ho
    %     plot(rrG)
    clear rater ntM ntG d dd
    save(di(dIter).name)
end
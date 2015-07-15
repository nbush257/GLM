% single variable correlations
clear
ca
c
D = dir('*GLM.mat')
binsize = 50;
for II = 1:length(D)
    II
    load(D(II).name)
    opts = statset('UseParallel',1);
    if iscolumn(spikevec); spikevec = spikevec';end
    if size(mech_85,1)<size(mech_85,2)
        mech_85 = mech_85';
        geo_85 = geo_85';
    end
    mech_85 = naninterp(mech_85);
    geo_85 = naninterp(geo_85);
    if exist('prox','var')
        C(prox) = 0;
    end
    mech_85(~C,:)=0;
    geo_85(~C,:)=0;
    
    rate_mech = tsmovavg(mech_85','s',binsize);rate_mech(:,isnan(rate_mech(1,:))) = [];
    rate_geo = tsmovavg(geo_85','s',binsize);rate_geo(:,isnan(rate_geo(1,:))) = [];
    rate_spike = tsmovavg(spikevec,'s',binsize);rate_spike(isnan(rate_spike)) = [];
    
    if II == 1
        keeperMech = rate_mech;
        keeperGeo = rate_geo;
        keeperRate = rate_spike;
    end
    
    mechCorr = bootstrp(100,@corr,rate_mech',rate_spike','Options',opts);
    mechCorr = mechCorr.^2;
    mechSE = std(mechCorr);
    geoCorr = bootstrp(100,@corr,rate_geo',rate_spike','Options',opts);
    geoCorr = geoCorr.^2;
    geoSE = std(geoCorr);
    Corr(II).FX = mechCorr(:,1);
    Corr(II).FY = mechCorr(:,2);
    Corr(II).M = mechCorr(:,3);
    
    Corr(II).R = geoCorr(:,1);
    Corr(II).TH = geoCorr(:,2);
    
    SE(II).FX = mechSE(:,1);
    SE(II).FY = mechSE(:,2);
    SE(II).M = mechSE(:,3);
    
    SE(II).R = geoSE(:,1);
    SE(II).TH = geoSE(:,2);
    
end
clearvars -except Corr SE
errorbar(mean([Corr.M]),[SE.M],'ko');ho
errorbar(mean([Corr.FX]),[SE.FX],'ko')
errorbar(mean([Corr.TH]),[SE.R],'ko')
errorbar(mean([Corr.R]),[SE.TH],'ko')


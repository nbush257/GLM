% single variable correlations
D = dir('*GLM.mat')
binsize = 100;
clear p RSQ
for II = 1:length(D)
    II
    load(D(II).name)
    if iscolumn(spikevec); spikevec = spikevec';end
    if size(mech_85,2)<size(mech_85,1)
        mech_85 = mech_85';
        geo_85 = geo_85';
    end
    mech_85 = naninterp(mech_85);
    geo_85 = naninterp(geo_85);
    mech_85(:,~C)=0;
    geo_85(:,~C)=0;
    
    rate_mech = tsmovavg(mech_85,'s',binsize);
    rate_geo = tsmovavg(geo_85,'s',binsize);
    rate_spike = tsmovavg(spikevec,'s',binsize);
    
    FXfit = fitlm(rate_mech(1,:),rate_spike);
    FYfit = fitlm(rate_mech(2,:),rate_spike);
    Mfit = fitlm(rate_mech(3,:),rate_spike);
    
    Rfit = fitlm(rate_geo(1,:),rate_spike);
    THfit = fitlm(rate_geo(2,:),rate_spike);
    
    p(II).FX = FXfit.coefTest;
    p(II).FY = FYfit.coefTest;
    p(II).M = Mfit.coefTest;
    p(II).R = Rfit.coefTest;
    p(II).TH = THfit.coefTest;
    
    RSQ(II).FX = FXfit.Rsquared.Adjusted;
    RSQ(II).FY = FYfit.Rsquared.Adjusted;
    RSQ(II).M = Mfit.Rsquared.Adjusted;
    RSQ(II).R = Rfit.Rsquared.Adjusted;
    RSQ(II).TH = THfit.Rsquared.Adjusted;
    
    
end
clearvars -except p RSQ
%% plot
rmR = [p.R]>(.05);
rmFX = [p.FX]>(.05);
rmFY = [p.FY]>(.05);
rmM = [p.M]>(.05);
rmTH = [p.TH]>(.05);

R = [RSQ.R];
FX = [RSQ.FX];
FY = [RSQ.FY];
M = [RSQ.M];
TH = [RSQ.TH];

RNS = R(rmR);
FXNS = FX(rmFX);
FYNS = FY(rmFY);
MNS = M(rmM);
THNS = TH(rmTH);

R(rmR) = NaN;
FX(rmFX) = NaN;
FY(rmFY) = NaN;
M(rmM) = NaN;
TH(rmTH) = NaN;

ca
fig
plot(R,'*');ho;plot(find(rmR),repmat(-.1,sum(rmR),1),'ko');title('R')
axy(-.2,1)
hline(0)


fig
plot(TH,'*');ho;plot(find(rmTH),repmat(-.1,sum(rmTH),1),'ko');title('TH')
axy(-.2,1)
hline(0)

fig
plot(FX,'*');ho;plot(find(rmFX),repmat(-.1,sum(rmFX),1),'ko');title('FX')
axy(-.2,1)
hline(0)

fig
plot(FY,'*');ho;plot(find(rmFY),repmat(-.1,sum(rmFY),1),'ko');title('FY')
axy(-.2,1)
hline(0)

fig
plot(M,'*');ho;plot(find(rmM),repmat(-.1,sum(rmM),1),'ko');title('M')
axy(-.2,1)
hline(0)







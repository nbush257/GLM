d = dir('*Hist.mat');
% Names
names = {}
for ii = 1:length(d)
    if ~isempty(regexp(d(ii).name,'Rat1*'))
        names{ii} = d(ii).name(9:16)
    else
        names{ii} = d(ii).name([9:10 21 22 28:33])
    end
end
% compare ROC
for ii = 1:37
    subplot(7,6,ii)
    plot(geo(ii).ROC_x,geo(ii).ROC_y)
    ho
    plot(mech(ii).ROC_x,mech(ii).ROC_y)
    ho
    title(names{ii},'Interpreter','None')
    dline
end
% compare AUC
subplot(221)
plot([[mech.AUC];[geo.AUC]],'k')
set(gca,'XTickLabel',{'Mech','Geo'})
title('All')
ylabel('AUC')
subplot(222)
plot([[m_No_Deriv.AUC];[g_No_Deriv.AUC]],'k')
title('No Derivative')
set(gca,'XTickLabel',{'Mech','Geo'})
ylabel('AUC')
subplot(223)
plot([[m_JustDeriv.AUC];[g_JustDeriv.AUC]],'k')
title('Just Derivative')
ylabel('AUC')
set(gca,'XTickLabel',{'Mech','Geo'})


%%
figure
mG = mean([geo.AUC]);sG = std([geo.AUC])/sqrt(numel([geo.AUC]));
mGJ = mean([g_JustDeriv.AUC]);sGJ = std([g_JustDeriv.AUC])/sqrt(numel([g_JustDeriv.AUC]));
mGN = mean([g_No_Deriv.AUC]);sGN = std([g_No_Deriv.AUC])/sqrt(numel([g_No_Deriv.AUC]));

mM = mean([mech.AUC]);sM = std([mech.AUC])/sqrt(numel([mech.AUC]));
mMJ = mean([m_JustDeriv.AUC]);sMJ = std([m_JustDeriv.AUC])/sqrt(numel([m_JustDeriv.AUC]));
mMN = mean([m_No_Deriv.AUC]);sMN = std([m_No_Deriv.AUC])/sqrt(numel([m_No_Deriv.AUC]));
%bar([mG mM;mGJ mMJ; mGN mMN]);
ho
errorbar_groups([mG mGJ mGN;mM mMJ mMN],[sG sGJ sGN ;sM sMJ sMN])

set(gca,'XTickLabel',{'All','Just Deriv','No Deriv'})
legend({'Geo','Mech'})
title('Area Under ROC')
hline(.5)
set(gca,'FontSize',20)

%% rates
rate = struct;
rate.mech = [];
rate.geo = [];
rate.spike = [];
bin = 5;

for ii = 1:37
    rate(ii).mech = tsmovavg(mech(ii).P','s',bin);
    rate(ii).geo= tsmovavg(geo(ii).P','s',bin);
    rate(ii).spike = tsmovavg(mech(ii).tested_spikes','s',bin);
    rate(ii).geo(isnan(rate(ii).geo)) = 0;
    rate(ii).spike(isnan(rate(ii).spike)) = 0;
    rate(ii).mech(isnan(rate(ii).mech)) = 0;
    
end

%% Predictive mode
figure
for ii = 1:length(geo)
    ii
    [G(ii).x,G(ii).y,~,G(ii).AUC] = perfcurve(geo(ii).tested_spikes,geo(ii).predictiveModeY,'1');
    [M(ii).x,M(ii).y,~,M(ii).AUC] = perfcurve(mech(ii).tested_spikes,mech(ii).predictiveModeY,'1');
    subplot(6,7,ii)
    plot(G(ii).x,G(ii).y)
    ho
    plot(M(ii).x,M(ii).y)
    title(names{ii})
end

%% No Hist
nl = @(x) 1./(1 +exp(-x));
for ii = 1:length(geo)
geoLogit{ii} = nl(geo(ii).Y);
mechLogit{ii} = nl(mech(ii).Y);
end

for ii = 1:length(geo)
    ii
    [G_noHist(ii).x,G_noHist(ii).y,~,G_noHist(ii).AUC] = perfcurve(geo(ii).tested_spikes,geoLogit{ii},'1');
    [M_noHist(ii).x,M_noHist(ii).y,~,M_noHist(ii).AUC] = perfcurve(mech(ii).tested_spikes,mechLogit{ii},'1');
end

figure
for ii = 1:length(geo)
    subplot(6,7,ii)
    plot(G_noHist(ii).x,G_noHist(ii).y)
    ho
    plot(M_noHist(ii).x,M_noHist(ii).y)
    title(names{ii})
end

mG_noHist = mean([G_noHist.AUC]);  sG_noHist = std([G_noHist.AUC])/sqrt(numel([G_noHist.AUC]));
mM_noHist = mean([M_noHist.AUC]);  sM_noHist = std([M_noHist.AUC])/sqrt(numel([M_noHist.AUC]));
errorbar([mG_noHist mM_noHist],[sG_noHist sM_noHist]);
plot([[G_noHist.AUC];[M_noHist.AUC]],'k')

%% med v dis


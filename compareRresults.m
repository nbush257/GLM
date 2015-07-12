% compare ROC
for ii = 1:37
    subplot(7,6,ii)
    plot(geo(ii).ROC_x,geo(ii).ROC_y)
    ho
    plot(mech(ii).ROC_x,mech(ii).ROC_y)
    ho
    dline
end
% compare AUC
figure
plot([[mech.AUC];[geo.AUC]],'k')
figure
plot([[m_No_Deriv.AUC];[g_No_Deriv.AUC]],'k')



figure
mG = mean([geo.AUC]);
mGJ = mean([g_JustDeriv.AUC]);
mGN = mean([g_No_Deriv.AUC]);

mM = mean([mech.AUC]);
mMJ = mean([m_JustDeriv.AUC]);
mMN = mean([m_No_Deriv.AUC]);
bar([mG mM;mGJ mMJ; mGN mMN])
set(gca,'XTickLabel',{'All','Just Deriv','No Deriv'})
legend({'Geo','Mech'})
title('Area Under ROC')

%% rates
rate = struct;
rate.mech = [];
rate.geo = [];
rate.spike = [];


for ii = 1:37
    rate(ii).mech = tsmovavg(mech(ii).P','s',15);
    rate(ii).geo= tsmovavg(geo(ii).P','s',15);
    rate(ii).spike = tsmovavg(mech(ii).tested_spikes','s',15);
    rate(ii).geo(isnan(rate(ii).geo)) = 0;
    rate(ii).spike(isnan(rate(ii).spike)) = 0;
    rate(ii).mech(isnan(rate(ii).mech)) = 0;
    
end

close all
%% Plot Awake versus anaesthetized
f1 = figure
plot([rG_all(1:10);rM_all(1:10)],'k-.');ho
plot([rG_all(11:end);rM_all(11:end)],'r-.')
mW = plot([mean(rG_all(1:10)) mean(rM_all(1:10))]','k-o','DisplayName','Awake');ln3
mA = plot([nanmean(rG_all(11:end)) nanmean(rM_all(11:end))]','r-*','DisplayName','Anaesthetized');ln3
ax  = gca
ax.XTick = [1 2]
ax.XTickLabel = {'Geo','Mech'}
axx(0.5,2.5)
axy(0,1)
legend([mW mA])
title('Awake Versus Anaesthetized')
ylabel('Pearson Correlation Coefficient')

%% Plot medial versus distal
f2 = figure;
plot([rG_med_all;rG_dis_all],'r-.');ho;
plot([rM_med_all;rM_dis_all],'k-.')
lG = plot([nanmean(rG_med_all) nanmean(rG_dis_all)]','r-o','MarkerFaceColor','r','DisplayName','Geo');ln3
lM = plot([nanmean(rM_med_all) nanmean(rM_dis_all)]','k-o','MarkerFaceColor','k','DisplayName','Mech');ln3
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel = {'Med','Dis'};
legend([lG lM])
axx(.5,2.5)
axy(0,1)
title('Medial Versus Distal')
ylabel('Pearson Correlation Coefficient')
%% plot PCC versus preferred direction
f3 = figure;
dis2Pref = abs(cos(deg2rad(allprefDir)))

plot(dis2Pref,rG_all,'r*')
ho
plot(dis2Pref,rM_all,'k*')
axx(-.2,1.2)
xlabel('abs(Cos(PD))')
ylabel('Pearson Correlation Coefficient')
ybG = polyfit(dis2Pref,rG_all,1);
ybM = polyfit(dis2Pref,rM_all,1);

yG = polyval(ybG,[0 1]);
yM = polyval(ybM,[0 1]);
pdG = plot([0 1],yG,'r-o','DisplayName','Geo');
pdM  = plot([0 1],yM,'k-o','DisplayName','Mech');
ax = gca;
ax.XTick=[0 1];
legend([pdG pdM]);
axy(0,1)
title('Preferred Direction')
%%
saveas(f1,'awakeVanaesthetized.fig','fig')
saveas(f2,'medVdis.fig','fig')
saveas(f3,'PD.fig','fig')

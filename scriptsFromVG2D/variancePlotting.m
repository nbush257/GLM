ca
Stim = [mech_85 geo_85 vel];
aa = [];
filters = weights.X.data;



% subplots
for ii = 1:6
    aa(:,ii) = conv(Stim(:,ii),filters(:,ii));
end
aa(aa==0)=NaN;

figure
ho
for ii = 1:6
    histogram(aa(:,ii))
end
legend(Names)
figure
for ii = 1:6
    
    ha(ii) = subplot(2,3,ii);
    histogram(aa(:,ii));
    title(Names{ii})
    
end


linkaxes(ha);
% legend(Names);
figure
opt = statset('UseParallel',true);
varBoot = bootstrp(100,@nanvar,aa,'Options',opt);
ho

bar(mean(varBoot)); errorbar(mean(varBoot),std(varBoot),'r.','MarkerSize',10,'LineStyle','none')
title('Variance of the Filtered Signal')
set(gca,'XTick',[1,2,3,4,5,6],'XtickLabel',Names)
ylabel('Variance')
figure
plot(aa)
legend(Names)
title('Filtered Signal Over Time')
xlabel('Time (ms)')
ylabel('Filtered signal value')


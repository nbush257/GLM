bnh = bnh;
filters = buildGLM.combineWeights(dm_nh,bnh(2:end));
filters = filters.X.data;
Stim = [mech_85 geo_85 vel];
aa = zeros(size(Stim));
for ii = 1:size(Stim,2)
    aa(:,ii) = conv(Stim(:,ii),filters(:,ii),'same');
end
aa(aa==0) = NaN;

ho
for ii = 1:6

    histogram(aa(:,ii))
    
end
legend(Names)
figure
bar(nanmedian(abs(aa)));ho
errorbar(1:6, nanmedian(abs(aa)),prctile(abs(aa),25),prctile(abs(aa),75),'r.','LineWidth',3)
set(gca,'XTick',[1 2 3 4 5 6])
Names2 = Names;
Names2{5} = '\theta';
Names2{1} = 'F_x';
Names2{2} = 'F_y';

set(gca,'XTickLabel',Names2)

dFilt = dir('*.mat')
fBoth = figure;
fMech = figure;
fGeo = figure;
colors = {'k','b','r','g','m'};


for dd = 1:length(dFilt)
    load(dFilt(dd).name,'weight*')
    
    if any(regexp(dFilt(dd).name,'_N\d{2}'))
        idx = regexp(dFilt(dd).name,'N\d{2}');
        namer = dFilt(dd).name(idx:idx+2);
        
    else
        namer = dFilt(dd).name([9:10 21:22 28:33]);
    end
    for ii = 1:5
        bWeight(ii,:,:) = weightB(ii).X.data;
        mWeight(ii,:,:) = weightM(ii).X.data;
        gWeight(ii,:,:) = weightG(ii).X.data;
    end
    mBweight = squeeze(mean(bWeight));
    seBweight = squeeze(std(bWeight))/sqrt(size(bWeight,1));
    
    mMweight = squeeze(mean(mWeight));
    seMweight = squeeze(std(mWeight))/sqrt(size(mWeight,1));
    
    mGweight = squeeze(mean(gWeight));
    seGweight = squeeze(std(gWeight))/sqrt(size(gWeight,1));
    %%
    figure(fBoth)
    
    subplot(6,6,dd)
    ho
    for ii =1:5
        shadedErrorBar(weightB(1).X.tr(:,ii),mBweight(:,ii),seBweight(:,ii),colors{ii})
    end
    title(namer)
    hold off
    %%
    figure(fMech)
    subplot(6,6,dd)
    title(namer)
    hold on   
    for ii = 1:3
        shadedErrorBar(weightM(1).X.tr(:,ii),mMweight(:,ii),seMweight(:,ii),colors{ii})
    end
    hold off
    %%
    figure(fGeo)
    subplot(6,6,dd)
    title(namer)
    hold on
    for ii = 1:2
        shadedErrorBar(weightG(1).X.tr(:,ii),mGweight(:,ii),seGweight(:,ii),colors{ii+3})
    end
    hold off
end
figure(fBoth)
subplot(6,6,dd+1)
for ii = 1:5
    ho
plot([1 2],[ii ii],colors{ii})
end
legend({'FX','FY','M','R','\Theta'})
title('Legend')
hold off

figure(fMech)
hold on
subplot(6,6,dd+1)
title('Legend')
for ii = 1:3
    ho
plot([1 2],[ii ii],colors{ii})
end
legend({'FX','FY','M'})
hold off
figure(fGeo)
subplot(6,6,dd+1)
for ii = 1:2
    ho
plot([1 2],[ii ii],colors{ii+3})
end
title('Legend')
legend({'R','\Theta'})
hold off


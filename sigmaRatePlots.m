names = dir('*noM.mat')
for ii = 1:length(names)
    if any(regexp(names(ii).name,'Rat1*'))
        nTitle{ii} = names(ii).name(9:11);
    else
        nTitle{ii} = names(ii).name([9:10 21 22 28:33])
    end
end
%% Find ISI peaks
idx = [];
for ii = 1:37
    ISI = diff(find(spikes{ii}));
    ISI(ISI>200) = [];
    [N,X] = hist(ISI,100);
    [mN,idmN] = max(N);
    idx(ii) = round(X(find(N<(max(N)/2) & X> X(idmN),1)));
    idx(ii)
    
end
%% manually grab the knee_pt
% for ii = 1:length(r)
%     ca
%     plot(stdRange(1:150),r(ii).noHistG(1:150));ln2
%     ho
%     plot(stdRange(1:150),r(ii).noHistM(1:150));ln2
%
%     [x,~] = ginput(1);
%     gIdx(ii) = dsearchn([r(ii).noHistG(1:150)]',[x]);
%     mIdx(ii) = dsearchn([r(ii).noHistM(1:150)]',[x]);
%     cla
%
% end

%% compare derivative
f1 = figure;
f2 = figure;
for ii= 1:37
    figure(f1)
    subplot(6,7,ii)
    ho
    plot(r(ii).noHistG(1:200),'g');ln2
    plot(r(ii).noHistGND(1:200),'r');ln2
    plot(r(ii).noHistGJD(1:200),'b');ln2
    title(nTitle{ii})
    plot(idx(ii),r(ii).noHistG(idx(ii)),'k*')
    plot(idx(ii),r(ii).noHistGND(idx(ii)),'k*')
    plot(idx(ii),r(ii).noHistGJD(idx(ii)),'k*')
    
    
    figure(f2)
    subplot(6,7,ii)
    ho
    plot(r(ii).noHistM(1:200),'g');ln2
    plot(r(ii).noHistMND(1:200),'r');ln2
    plot(r(ii).noHistMJD(1:200),'b');ln2
    
    plot(idx(ii),r(ii).noHistM(idx(ii)),'k*')
    plot(idx(ii),r(ii).noHistMND(idx(ii)),'k*')
    plot(idx(ii),r(ii).noHistMJD(idx(ii)),'k*')
    title(nTitle{ii})
    
end
figure(f1)
subplot(6,7,ii+1)
plot([0 1],[1 1],'g');ho
plot([0 1],[2 2],'r')
plot([0 1],[3 3],'b')
legend({'All','No Deriv','Just Deriv'})
title('Geometry')

figure(f2)
subplot(6,7,ii+1)
plot([0 1],[1 1],'g');ho
plot([0 1],[2 2],'r')
plot([0 1],[3 3],'b')
legend({'All','No Deriv','Just Deriv'})
title('mechanics')
%% compare model types

f1 = figure;
f2 = figure;
n = fieldnames(r);
for ii = 1:37
    figure(f1)
    subplot(6,7,ii)
    plot(r(ii).generativeG(1:200),'b');ln2
    ho
    plot(r(ii).predictiveG(1:200),'r');ln2
    plot(r(ii).noHistG(1:200),'g');ln2
    
    plot(idx(ii),r(ii).generativeG(idx(ii)),'k*')
    plot(idx(ii),r(ii).predictiveG(idx(ii)),'k*')
    plot(idx(ii),r(ii).noHistG(idx(ii)),'k*')
    title(nTitle{ii})
    
    figure(f2)
    subplot(6,7,ii)
    plot(r(ii).generativeM(1:200),'b');ln2
    ho
    plot(r(ii).predictiveM(1:200),'r');ln2
    plot(r(ii).noHistM(1:200),'g');ln2
    
    plot(idx(ii),r(ii).generativeM(idx(ii)),'k*')
    plot(idx(ii),r(ii).predictiveM(idx(ii)),'k*')
    plot(idx(ii),r(ii).noHistM(idx(ii)),'k*')
    
    title(nTitle{ii})
end

figure(f1)
subplot(6,7,ii+1)
plot([0 1],[1 1],'b');ho
plot([0 1],[2 2],'r')
plot([0 1],[3 3],'g')
legend({'generative','predictive','noHist'})
title('Geometry')

figure(f2)
subplot(6,7,ii+1)
plot([0 1],[1 1],'b');ho
plot([0 1],[2 2],'r')
plot([0 1],[3 3],'g')
legend({'generative','predictive','noHist'})
title('mechanics')
%% compare geo and mech
toCompare = {'noHistG','noHistM'};
color = {'g','r'};
figure
for ii =1:37
    subplot(6,7,ii)
    
    for jj = 1:length(toCompare)
        plot(r(ii).(toCompare{jj}),color{jj});ln2;ho
        plot(idx(ii),r(ii).(toCompare{jj})(idx(ii)),'k*')
        title(nTitle{ii})
        
    end
    for jj = 1:length(toCompare)
        plot(idx(ii),r(ii).(toCompare{jj})(idx(ii)),'k*')
    end
    
    
end

subplot(6,7,ii+1)
for ii = 1:length(toCompare)
    plot([0 1],[ii ii],color{ii});ho;
end
legend(toCompare)

%%

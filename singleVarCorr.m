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
    C(1) = 0;
    C(end) = 0;
    mech_85 = naninterp(mech_85);
    geo_85 = naninterp(geo_85);
    if exist('prox','var')
        C(prox) = 0;
    end
    mech_85(~C,:)=0;
    geo_85(~C,:)=0;
    
    %     rate_mech = tsmovavg(mech_85','s',binsize);rate_mech(:,isnan(rate_mech(1,:))) = [];
    %     rate_geo = tsmovavg(geo_85','s',binsize);rate_geo(:,isnan(rate_geo(1,:))) = [];
    rate_spike = smoothts(spikevec,'g',length(spikevec),15);%tsmovavg(spikevec,'s',binsize);rate_spike(isnan(rate_spike)) = [];
    
    starts = find(diff(C)==1)+1;
    stops = find(diff(C)==-1);
    newMech = [];newGeo = [];newSpike= [];
    for ii = 1:length(starts)
        newMech = [newMech;zeros(50,3);mech_85(starts(ii):stops(ii),:);zeros(50,3)];
        newGeo = [newGeo;zeros(50,2);geo_85(starts(ii):stops(ii),:);zeros(50,2)];
        newSpike = [newSpike zeros(1,50) rate_spike(starts(ii):stops(ii)) zeros(1,50)];
    end
    mech_85 = newMech;
    geo_85 = newGeo;
    rate_spike = newSpike;
    vel = geo_85(:,2);
    vel = cdiff(vel)';
    vel(isnan(vel))=0;
    
    if II == 1
        keeperMech = mech_85;
        keeperGeo = [geo_85 vel];
        keeperRate = rate_spike;
        
    end
    tic
    outCorr = bootstrp(500,@corr,[mech_85 geo_85 vel],rate_spike','Options',opts);
    b = toc;
    fprintf('Correlations took %.4f seconds for %i ms\n',b,length(spikevec))
    Corr(II).FX = outCorr(:,1);
    Corr(II).FY = outCorr(:,2);
    Corr(II).M = outCorr(:,3);
    
    Corr(II).R = outCorr(:,4);
    Corr(II).TH = outCorr(:,5);
    Corr(II).V = outCorr(:,6);
    %
    %     CI(II).FX = outCI(:,1);
    %     CI(II).FY = outCI(:,2);
    %     CI(II).M = outCI(:,3);
    %
    %     CI(II).R = outCI(:,4);
    %     CI(II).TH = outCI(:,5);
    clearvars -except Corr CI keeper* D II binsize
    tic
    %save('Corr.mat','Corr','keeper*')
    b = toc;
    if b>1
        fprintf('Saving took %f seconds\n',b)
    end
end
clearvars -except Corr keeper*
names = fieldnames(Corr);
for ii = 1:length(Corr)
    for jj = 1:length(names)
        erUB(ii).(names{jj}) = prctile(Corr(ii).(names{jj}),97.5) - mean(Corr(ii).(names{jj}));
        erLB(ii).(names{jj}) = -prctile(Corr(ii).(names{jj}),2.5) + mean(Corr(ii).(names{jj}));
    end
end
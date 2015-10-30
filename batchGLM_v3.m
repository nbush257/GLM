clear
ca
d = dir('*toGLM.mat');
noProx = 1;
saveTGL =1;% 1 if you want to save the outputs
basisSize = 2;
lastPeak = 5;
nlOffset = 4;
nHist = 3;
k=5;

suffix = 'scaleGLM';
saveLoc = 'postRev';mkdir(saveLoc)
stimBasis = basisFactory.makeNonlinearRaisedCos(basisSize,3,[0 lastPeak],nlOffset);
plot(stimBasis.B)
clearvars -except d ddd cellType noProx k nHist lastPeak basisSize nlOffset rateBin histSize saveTGL suffix stimBasis histBasis saveLoc
d = dir('*toGLM.mat');

for ddd = 1:length(d)
    fname = [d(ddd).name(1:end-10) suffix];% output filename -no suffix
    clear med prox dis
    load(d(ddd).name,'mech*','geo*','C','spike*','med','prox','dis')
    nl = @(x) 1./(1 +exp(-x));
    C = logical(C);
    C(1) = 0;
    C(end) = 0;
    
    if size(mech_85,2)>size(mech_85,1)
        mech_85 = mech_85';
        geo_85 = geo_85';
    end
    
    if sum(spikevec)<50
        continue
    end
    if isrow(spikevec);spikevec = spikevec'; end;
    if isrow(C);C = C'; end;
    
    if noProx & exist('prox','var')
        C(prox) = [];
        mech_85(prox,:) = [];
        geo_85(prox,:) = [];
        spikevec(prox) = [];
        med(prox) = [];
        dis(prox) = [];
    end
    cstart = find(diff(C)==1)+1;
    cend = find(diff(C)==-1);
    for iii = 1:length(cstart)
        
        tM = mech_85(cstart(iii):cend(iii),:);
        
        fM = tM(~isnan(tM(:,1)),:);
        %         fG = tG(~isnan(tG(:,1)),:);
        if ~isempty(fM)
            fM = fM(1,:);
        else
            fM = [NaN NaN NaN];
        end
        %         if ~isempty(fG)
        %             fG = fG(1,:);
        %         else
        %             fG = [NaN NaN];
        %         end
        
        
        
        %         geo_85(cstart(iii):cend(iii),:) = geo_85(cstart(iii):cend(iii),:)-repmat(fG,cend(iii)-cstart(iii)+1,1);
        mech_85(cstart(iii):cend(iii),:) = mech_85(cstart(iii):cend(iii),:)-repmat(fM,cend(iii)-cstart(iii)+1,1);
    end
    
    
    
    vel = cdiff(geo_85(:,2))';
    
    %% scaling
    fprintf('scaling...\n')
    mech_85 = mech_85./repmat(nanstd(mech_85),length(spikevec),1);
    geo_85= geo_85./repmat(nanstd(geo_85),length(spikevec),1);
    vel = vel./nanstd(vel);
    
    %%
    
    replacerM = isnan(mech_85);
    replacerG = isnan(geo_85);
    replacerVel = isnan(vel);
    mech_85(isnan(mech_85)) = 0;
    geo_85(isnan(geo_85)) = 0;
    vel(isnan(vel)) = 0;
    
    mech_85(replacerM) = 0;
    
    geo_85(replacerG) = 0;
    
    vel(replacerVel) = 0;
    
    %     %%
    %     fprintf('doubling into pos and neg')
    %     mPos = mech_85; mNeg = mech_85;
    %     mPos(mPos<0) = 0; mNeg(mNeg>0) = 0;
    %     gPos = geo_85; gNeg = geo_85;
    %     gPos(gPos<0) = 0; gNeg(gNeg>0) = 0;
    
    
    
    
    [XM_nh,dm_nh] = buildDesignMatrix([mech_85 geo_85 vel],spikevec,'deriv',0,'hist',0,'bsStim',stimBasis);
    XM_h = XM_nh;
    for iii = 1:nHist
        XM_h = [XM_h [zeros(iii,1);spikevec(1:end-iii)]];
    end
    yHatH = zeros(size(spikevec));
    yHatSimIn = zeros(size(spikevec));
    %     cvnet = cvglmnet(XM_nh,spikevec,'binomial');
    %     idx = crossvalind('Kfold',length(spikevec),k);
    %     for kk = 1:k
    %         [H(ddd,kk).b,~,H(ddd,kk).stats] = glmfit(XM_h(idx~=kk,:),spikevec(idx~=kk),'binomial');
    %         [NH(ddd,kk).b,~,NH(ddd,kk).stats] = glmfit(XM_nh(idx~=kk,:),spikevec(idx~=kk),'binomial');
    %         yHatSimIn(idx==kk) = glmval(H(ddd,kk).b(1:end-nHist),XM_nh(idx==kk,:),'identity');
    %         yHatH(idx==kk) = glmval(H(ddd,kk).b,XM_h(idx==kk,:),'logit');
    %         yHatNH(idx==kk) = glmval(NH(ddd,kk).b,XM_nh(idx==kk,:),'logit');
    %
    % %     end
    %     fprintf('Simulating responses')
    %     hist = mean([H(ddd,:).b],2);hist = hist(end-nHist+1:end);
    %     ySim = simGLM4(yHatSimIn,hist,10000);
    % [mB{ddd},mDev{ddd},mStats{ddd}] = glmfit(XM_h(:,[1:6 end-nHist+1:end]),spikevec,'binomial');
    % [gB{ddd},gDev{ddd},gStats{ddd}] = glmfit(XM_h(:,7:end),spikevec,'binomial');
    % [B{ddd},dev{ddd},stats{ddd}] = glmfit(XM_h,spikevec,'binomial');
    %
    %
    % [mBnh{ddd},mDevnh{ddd},mStatsnh{ddd}] = glmfit(XM_nh(:,[1:6]),spikevec,'binomial');
    % [gBnh{ddd},gDevnh{ddd},gStatsnh{ddd}] = glmfit(XM_nh(:,7:end),spikevec,'binomial');
    % [Bnh{ddd},devnh{ddd},statsnh{ddd}] = glmfit(XM_nh,spikevec,'binomial');
    %
    % testSpikes{ddd} = spikevec(C);
    % testXM_nh{ddd} = XM_nh(C,:);
    % testXM_h{ddd} = XM_h(C,:);
    %
    [bnh,~,stats_nh] = glmfit(XM_nh,spikevec,'binomial');
    [bh,~,stats_h] = glmfit(XM_h,spikevec,'binomial');
    yHat = glmval(bh(1:end-nHist),XM_nh,'identity');
    fprintf('simulating...\n')
    ySim = simGLM4(yHat,bh(end-nHist+1:end),1000);
    yHat_nh = glmval(bnh,XM_nh,'logit');
    yHat_hist = glmval(bh,XM_h,'logit');
    %
    % B{ddd} = bh;
    % BNH{ddd} = bnh;
    % DMNH{ddd} = dm_nh;
    % WEIGHTS{ddd} = buildGLM.combineWeights(dm_nh,bh(2:end-nHist));
    Names = {'FX','FY','M','R','THETA','V'}
    % YSIM{ddd} = ySim;
    % YHAT_NOHIST{ddd} = glmval(bnh,XM_nh,'logit');
    % allXM_H{ddd} = XM_h;
    % allXM_NH{ddd} = XM_nh;
    if ~exist('prox')
        prox = [];
        med = [];
        dis = [];
    end
    if saveTGL
        cd(saveLoc)
        save([d(ddd).name([1:end-9]) suffix '.mat'],'bnh','stats_nh','bh','stats_h','yHat_nh','ySim','yHat_hist','Names','stimBasis','spikevec','C','prox','med','dis','dm_nh')
        cd ..
    end
end

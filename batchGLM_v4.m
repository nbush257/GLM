% version 4 is designed for use with the new median filtered data
clear
ca
d = dir('*w_phase.mat');
noProx = 1;
noHolds = 0;
saveTGL =1;% 1 if you want to save the outputs
basisSize = 5;
lastPeak =50/3;
nlOffset = 1;
doubling = 0;% this is not working in this version
histTGL = 0;
nHist = 0;
k=10;
filt = 150;

suffix = 'done_150Filt';
saveLoc = 'done';mkdir(saveLoc)
stimBasis = basisFactory.makeNonlinearRaisedCos(basisSize,1,[0 lastPeak],nlOffset);

plot(stimBasis.B);drawnow

for ddd =1:length(d)
    fprintf('working on cell %i\n',ddd)
    fname = [d(ddd).name(1:end-10) suffix];% output filename -no suffix
    clear med prox dis spikevec Mech Geo C
    load(d(ddd).name)
    %     prox([111597:122146 192093:205831]) = 1;
    
    nl = @(x) 1./(1 +exp(-x));
    if exist('C','var') & length(C)~=size(Mech.filtAll,1)
        clear C
    end
    if ~exist('C','var')% if we don't have a C variable, estimate it by taking when all the mechanical variables are NaNs
        C = all(~isnan(Mech.filtAll'));
    end
    
    if isrow(spikevec);spikevec = spikevec'; end;
    if isrow(C);C = C'; end;
    
    if filt~=0
        for ii = 1:3
            Mech.filtAll(:,ii) = interpNaNFilt(Mech.filtAll(:,ii),1000,filt);
            Geo.filtAll(:,ii) = interpNaNFilt(Geo.filtAll(:,ii),1000,filt);
            
        end
    end
    
    % Edit C so that we have matching contact onsets and offsets
    C(isnan(C)) = 0;
    C = logical(C);
    C(1) = 0;
    C(end) = 0;
    % remove spikes outside of contact
    spikevec(~C) = 0;
    
    %     % merge the holds binary with the prox binary into just prox.
    %     if isrow(removeHolds);removeHolds = removeHolds';end
    %     if noProx & noHolds & exist('prox')
    %         if size(prox) == size(removeHolds)
    %         prox = removeHolds | prox;
    %         else
    %             clear prox med dis
    %         end
    %     end
    
    if exist('prox','var') && numel(prox)~=numel(spikevec)
        warning('Is This Still Happeining?')
        clear prox med dis
    end
    
    % if we have less than 50 spikes, do not analyze
    
    
    % remove the proximal contacts
    if noProx & exist('prox','var')
        C(prox) = [];
        Mech.filtAll(prox,:) = [];
        Geo.filtAll(prox,:) = [];
        spikevec(prox) = [];
        med(prox) = [];
        dis(prox) = [];
    end
    if sum(spikevec)<50
        continue
    end
    
    % set NaN to zero
    Mech.filtAll(isnan(Mech.filtAll)) = 0;
    Geo.filtAll(isnan(Geo.filtAll)) = 0;
    
    
    % maybe do this to apply a nonlinearity to direction
    if doubling
        fprintf('doubling into pos and neg')
        mPos = mech_85; mNeg = mech_85;
        mPos(mPos<0) = 0; mNeg(mNeg>0) = 0;
        gPos = geo_85; gNeg = geo_85;
        gPos(gPos<0) = 0; gNeg(gNeg>0) = 0;
        velPos = vel; velNeg = vel;
        velPos(velPos<0) = 0;velNeg(velNeg>0) = 0;
        
        mech_85 = [mPos mNeg];
        geo_85 = [gPos gNeg];
        vel = [velPos velNeg];
    end
    
    
    
    
    %% Build the design matrices
    [XM_Full,dm_Full] = buildDesignMatrix([Mech.filtAll Geo.filtAll],spikevec,'deriv',0,'hist',0,'bsStim',stimBasis);
    [XM_Mech,dm_Mech] = buildDesignMatrix([Mech.filtAll],spikevec,'deriv',0,'hist',0,'bsStim',stimBasis);
    [XM_Geo,dm_Geo] = buildDesignMatrix([Geo.filtAll],spikevec,'deriv',0,'hist',0,'bsStim',stimBasis);
    
    
    
    
    % incorporate history
    if histTGL
        
        for iii = 1:nHist
            XM_Full = [XM_Full [zeros(iii,1);spikevec(1:end-iii)]];
            XM_Mech = [XM_Mech [zeros(iii,1);spikevec(1:end-iii)]];
            XM_Geo = [XM_Geo [zeros(iii,1);spikevec(1:end-iii)]];
        end
    end
    %     yHatH = zeros(size(spikevec));
    %     yHatSimIn = zeros(size(spikevec));
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
    %% Perform fitting using crossvalidation.
    idx = crossvalind('Kfold',length(spikevec),k);
    yHat_Full = zeros(size(spikevec));
    yHat_Mech = zeros(size(spikevec));
    yHat_Geo = zeros(size(spikevec));
    for cross = 1:k
        [bFull,~,stats_Full] = glmfit(XM_Full(idx~=cross,:),spikevec(idx~=cross),'binomial');
        [bMech,~,stats_Mech] = glmfit(XM_Mech(idx~=cross,:),spikevec(idx~=cross),'binomial');
        [bGeo,~,stats_Geo] = glmfit(XM_Geo(idx~=cross,:),spikevec(idx~=cross),'binomial');
        
        yHat_Full(idx==cross)= glmval(bFull,XM_Full(idx==cross ,:),'logit');
        yHat_Mech(idx==cross)= glmval(bMech,XM_Mech(idx==cross,:),'logit');
        yHat_Geo(idx==cross)= glmval(bGeo,XM_Geo(idx==cross ,:),'logit');
        
        BFull(:,cross) = bFull;
        BMech(:,cross) = bMech;
        BGeo(:,cross) = bGeo;
    end
    
    
    %% simulate if you have a history term
    %     [bh,~,stats_h] = glmfit(XM_h,spikevec,'binomial');
    %     yHat = glmval(bh(1:end-nHist),XM_nh,'identity');
    %     fprintf('simulating...\n')
    %     ySim = simGLM4(yHat,bh(end-nHist+1:end),1000);
    %     yHat_nh = glmval(bnh,XM_nh,'logit');
    %     yHat_hist = glmval(bh,XM_h,'logit');
    
    %
    
    %% save outputs for every cell in one variable
    rate = smoothts(spikevec','g',length(spikevec),15);
    r(ddd,1) = corr(rate(C)',yHat_Full(C));
    r(ddd,2) = corr(rate(C)',yHat_Mech(C));
    r(ddd,3) = corr(rate(C)',yHat_Geo(C));
    Names = {'FX','FY','M','R','THETA','V'};
    
    if ~exist('prox')
        prox = [];
        med = [];
        dis = [];
    end
    
    if saveTGL
        cd(saveLoc)
        save([d(ddd).name([1:end-9]) suffix '.mat'],'C','prox','med','dis','Geo','Mech','spikevec','yHat*','XM*','dm*','stimBasis','filt','r')
        cd ..
    end
    
    
end

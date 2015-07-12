dter = dir('*toGLM.mat')
R = {}; TH_o = {}; FY = {};FX = {};M = {};FR={};
for di = 1:length(dter)
    di
    load(dter(di).name)
    if exist('prox','var')
        C(prox) = 0;
    end
    C(1) = 0;
    C(end) = 0;
    starts = find(diff(C)==1)+1;
    stops = find(diff(C)==-1);
    
    avgR = NaN(size(starts));rangeFX = NaN(size(starts)); rangeTH = NaN(size(starts)); rangeFY = NaN(size(starts));rangeM = NaN(size(starts));avgFR = zeros(size(starts));
    if size(geo_85,2)>size(geo_85,1)
        geo_85 = geo_85';
        mech_85 = mech_85';
    end
    
    for jj = 1:length(starts)
        avgR(jj) = nanmean(geo_85(starts(jj):stops(jj),1));
        [~,idxTH] = nanmax(abs(geo_85(starts(jj):stops(jj),2)));
        rangeTH(jj) = geo_85(starts(jj)+idxTH-1,2);
        
        [~,idxFX] = nanmax(abs(mech_85(starts(jj):stops(jj),1)));
        rangeFX(jj) = mech_85(starts(jj)+idxFX-1,1);
        
        [~,idxFY] = nanmax(abs(mech_85(starts(jj):stops(jj),2)));
        rangeFY(jj) = mech_85(starts(jj)+idxFY-1,2);
        
        [~,idxM] = nanmax(abs(mech_85(starts(jj):stops(jj),3)));
        rangeM(jj) = mech_85(starts(jj)+idxM -1,3);
        avgFR(jj) = nanmean(spikevec(starts(jj):stops(jj)));
    end
    if length(starts)>5
        rMdl = fitlm(avgR,avgFR);
        thMdl = fitlm(rangeTH,avgFR);
        fxMdl = fitlm(rangeFX,avgFR);
        fyMdl = fitlm(rangeFY,avgFR);
        mMdl = fitlm(rangeM,avgFR);
        
        Rsquared(di).R = rMdl.Rsquared.Adjusted;
        Rsquared(di).TH = thMdl.Rsquared.Adjusted;
        Rsquared(di).FX = fxMdl.Rsquared.Adjusted;
        Rsquared(di).FY = fyMdl.Rsquared.Adjusted;
        Rsquared(di).M = mMdl.Rsquared.Adjusted;
        
        p(di).R = rMdl.coefTest;
        p(di).TH = thMdl.coefTest;
        p(di).FX = fxMdl.coefTest;
        p(di).FY = fyMdl.coefTest;
        p(di).M = mMdl.coefTest;
    end
    
        R{di} = avgR;
        TH_o{di} = rangeTH;
        FX{di} = rangeFX;
        FY{di} = rangeFY;
        M{di} = rangeM;
        FR{di} = avgFR;
end
%% show that R is not valid in awake, but yes in anaesth
awakeP = [p(1:10).R]<(.005) % 1 of 10 bonferroni p<.05
anaesthP = [p(11:end).R]<(.05/30) % 11 of 30 bonferroni p<.05

awakeR = [Rsquared(1:10).R];
anaesthR = [Rsquared(11:end).R];

nanmean(awakeR(awakeP)) % one example, .484
nanmean(anaesthR(anaesthP))% 11 examples min= .05, max = .44 mean = .2 
% R alone is not an accurate predictor of Firing Rate in active or passive
% case

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
awakePR = [p(1:10).R]<(.005) % 1 of 10 bonferroni p<.05
anaesthPR = [p(11:end).R]<(.05/30) % 11 of 30 bonferroni p<.05

awakeRR = [Rsquared(1:10).R];
anaesthRR = [Rsquared(11:end).R];

nanmean(awakeRR(awakePR)) % one example, .484
nanmean(anaesthRR(anaesthPR))% 11 examples min= .05, max = .44 mean = .2 
% R alone is not an accurate predictor of Firing Rate in active or passive
% case
%% M
awakePM = [p(1:10).M]<.005
anaesthPM = [p(11:end).M]<(.05/30) 


awakeRM = [Rsquared(1:10).M]
anaesthRM = [Rsquared(11:end).M]

nanmean(awakeRM(awakePM))
nanmean(anaesthRM(anaesthPM))
%% TH


awakePTH= [p(1:10).TH]<.005
anaesthPTH = [p(11:end).TH]<(.05/30) 


awakeRTH = [Rsquared(1:10).TH]
anaesthRTH = [Rsquared(11:end).TH]

nanmean(awakeRTH(awakePTH))
nanmean(anaesthRTH(anaesthPTH))
%% FX


awakePFX= [p(1:10).FX]<.005
anaesthPFX = [p(11:end).FX]<(.05/30) 

awakeRFX = [Rsquared(1:10).FX]
anaesthRFX = [Rsquared(11:end).FX]

nanmean(awakeRFX(awakePFX))
nanmean(anaesthRFX(anaesthPFX))
%% FX


awakePFY= [p(1:10).FY]<.005
anaesthPFY = [p(11:end).FY]<(.05/30) 

awakeRFY = [Rsquared(1:10).FY]
anaesthRFY = [Rsquared(11:end).FY]

nanmean(awakeRFY(awakePFY))
nanmean(anaesthRFY(anaesthPFY))
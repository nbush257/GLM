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
    
    avgR = zeros(size(starts));rangeFX = zeros(size(starts)); rangeTH = zeros(size(starts)); rangeFY = zeros(size(starts));rangeM = zeros(size(starts));avgFR = zeros(size(starts));
    if size(geo_85,2)>size(geo_85,1)
        geo_85 = geo_85';
        mech_85 = mech_85';
    end
    
    for jj = 1:length(starts)
        avgR(jj) = nanmean(geo_85(starts(jj):stops(jj),1));
        rangeTH(jj) = nanmax(geo_85(starts(jj):stops(jj),2));
        rangeFX(jj) = nanmax(mech_85(starts(jj):stops(jj),1));
        rangeFY(jj) = nanmax(mech_85(starts(jj):stops(jj),2));
        rangeM(jj) = nanmax(mech_85(starts(jj):stops(jj),3));
        avgFR(jj) = nanmean(spikevec(starts(jj):stops(jj)));
    end
    R{di} = avgR;
    TH_o{di} = rangeTH;
    FX{di} = rangeFX;
    FY{di} = rangeFY;
    M{di} = rangeM;
    FR{di} = avgFR;
end
for ii = 1:40
    plot(zscore(R{ii}),zscore(FR{ii}),'.')
    ho
end

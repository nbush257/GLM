clear
ca
d = dir('*toGLM.mat');
noProx = 1;
saveTGL =1;% 1 if you want to save the outputs
basisSize = 2;
lastPeak = 6;
nlOffset = 4;

suffix = '_simGLM';
saveLoc = 'postRev';mkdir(saveLoc)
stimBasis = basisFactory.makeNonlinearRaisedCos(basisSize,3,[0 lastPeak],nlOffset);




for ddd = 25%1:length(d)
    clearvars -except d ddd cellType lastPeak basisSize nlOffset rateBin histSize saveTGL suffix stimBasis histBasis saveLoc
    d = dir('*toGLM.mat');
    fname = [d(ddd).name(1:end-10) suffix];% output filename -no suffix
    load(d(ddd).name,'mech*','geo*','C','spike*','med','prox','dis')
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
    
    [XM_nh,dm_nh] = buildDesignMatrix(mech_85,spikevec,'deriv',1,'hist',0,'bsStim',stimBasis,'bsSpike',histBasis);
    XM
    mD.X  = XM;
    mD.X_noHist = XM_nh;
    mD.dm  = dm;
    mD.dm_noHist  = dm_nh;
    mD.spikes = spikevec;
    mD.C = C;
    mD.filtSize = length(stimBasis.tr);
    mD.histSize = histBasis.edim;
    out_mD = procVG_GLM(mD);
    
    [XM,dm] = buildDesignMatrix(mech_85,spikevec,'deriv',0,'hist',1,'bsStim',stimBasis,'bsSpike',histBasis);
    [XM_nh,dm_nh] = buildDesignMatrix(mech_85,spikevec,'deriv',0,'hist',0,'bsStim',stimBasis,'bsSpike',histBasis);
    mND.X  = XM;
    mND.X_noHist = XM_nh;
    mND.dm  = dm;
    mND.dm_noHist  = dm_nh;
    mND.spikes = spikevec;
    mND.C = C;
    mND.filtSize = length(stimBasis.tr);
    mND.histSize = histBasis.edim;
    out_mND = procVG_GLM(mND);
    
    [XM,dm] = buildDesignMatrix(mech_85,spikevec,'deriv',1,'covar',0,'hist',1,'bsStim',stimBasis,'bsSpike',histBasis);
    [XM_nh,dm_nh] = buildDesignMatrix(mech_85,spikevec,'deriv',1,'covar',0,'hist',0,'bsStim',stimBasis,'bsSpike',histBasis);
    mJD.X  = XM;
    mJD.X_noHist = XM_nh;
    mJD.dm  = dm;
    mJD.dm_noHist  = dm_nh;
    mJD.spikes = spikevec;
    mJD.C = C;
    mJD.filtSize = length(stimBasis.tr);
    mJD.histSize = histBasis.edim;
    out_mJD = procVG_GLM(mJD);
    
    
    [XG,dg] = buildDesignMatrix(geo_85,spikevec,'deriv',1,'hist',1,'bsStim',stimBasis,'bsSpike',histBasis);
    [XG_nh,dg_nh] = buildDesignMatrix(geo_85,spikevec,'deriv',1,'hist',0,'bsStim',stimBasis,'bsSpike',histBasis);
    gD.X  = XG;
    gD.X_noHist = XG_nh;
    gD.dm  = dg;
    gD.dm_noHist  = dg_nh;
    gD.spikes = spikevec;
    gD.C = C;
    gD.filtSize = length(stimBasis.tr);
    gD.histSize = histBasis.edim;
    out_gD = procVG_GLM(gD);
    
    [XG,dg] = buildDesignMatrix(geo_85,spikevec,'deriv',0,'hist',1,'bsStim',stimBasis,'bsSpike',histBasis);
    [XG_nh,dg_nh] = buildDesignMatrix(geo_85,spikevec,'deriv',0,'hist',0,'bsStim',stimBasis,'bsSpike',histBasis);
    gND.X = XG;
    gND.X_noHist = XG_nh;
    gND.dm = dg;
    gND.dm_noHist  = dg_nh;
    gND.spikes = spikevec;
    gND.C = C;
    gND.filtSize = length(stimBasis.tr);
    gND.histSize = histBasis.edim;
    out_gND = procVG_GLM(gND);
    
    
    [XG,dg] = buildDesignMatrix(geo_85,spikevec,'deriv',1,'covar',0,'hist',1,'bsStim',stimBasis,'bsSpike',histBasis);
    [XG_nh,dg_nh] = buildDesignMatrix(geo_85,spikevec,'deriv',1,'covar',0,'hist',0,'bsStim',stimBasis,'bsSpike',histBasis);
    gJD.X  = XG;
    gJD.X_noHist = XG_nh;
    gJD.dm = dg;
    gJD.dm_noHist = dg_nh;
    gJD.spikes = spikevec;
    gJD.C = C;
    gJD.filtSize = length(stimBasis.tr);
    gJD.histSize = histBasis.edim;
    out_gJD = procVG_GLM(gJD);
    
    
    
    
    if saveTGL
        
        cd(saveLoc)
        if exist('prox','var')
            save([fname suffix '.mat'],'out*','stimBasis','histBasis','mD','mND','gD','gND','med*','prox*','dis*')
        else
            save([fname suffix '.mat'],'out*','stimBasis','histBasis','mD','mND','gD','gND')
        end
        cd ..
    end
    
    
    
    
    
    
end


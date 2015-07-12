% get all vars
d = dir('*Hist.mat')
for dd = 1:length(d)
    load(d(dd).name)
    dd
    
    
    out_gD.raster = rmfield(out_gD,'raster');
    out_gJD.raster = rmfield(out_gJD,'raster');
    out_gND.raster = rmfield(out_gND,'raster');
    out_mD.raster = rmfield(out_mD,'raster');
    out_mJD.raster = rmfield(out_mJD,'raster');
    out_mND.raster = rmfield(out_mND,'raster');
    
    
    
    geo(dd) = out_gD;
    mech(dd) = out_mD;
    g_No_Deriv(dd) = out_gND;
    m_No_Deriv(dd) = out_mND;
    g_JustDeriv(dd) = out_gJD;
    m_JustDeriv(dd) = out_mJD;
    if exist('med','var')
        allMed{dd} = med;
        allProx{dd} = prox;
        allDis{dd} = dis;
    end
end
save('SummaryIndep.mat','geo','mech','g_No_Deriv','m_No_Deriv','g_JustDeriv','m_JustDeriv','allMed','allProx','allDis','-v7.3')
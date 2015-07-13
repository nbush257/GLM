d = dir('*GLM.mat')
nl = @(x) 1./(1 +exp(-x));
binsize =25;
r_geo = nan(1,length(d));
r_mech = nan(1,length(d));
for ii = 1:length(d)
    ii
    load(d(ii).name,'out*')
    rate = tsmovavg(out_gND.tested_spikes','s',binsize);rate(isnan(rate))=0;
    
    m_geo = mean(out_gND.raster,2);
    m_geo = tsmovavg(m_geo','s',binsize);
    m_geo(isnan(m_geo)) = 0;
    
    m_mech = mean(out_mND.raster,2);
    m_mech = tsmovavg(m_mech','s',binsize);
    m_mech(isnan(m_mech)) = 0;
    
    r_geo(ii) = corr(m_geo',rate');
    r_mech(ii) = corr(m_mech',rate');
    
    rate_geo{ii} = m_geo;
    rate_mech{ii} = m_mech;
    rate_spike{ii} = rate;
end
save('rateAndCorr.mat','rate_*','r_*','binsize')
    
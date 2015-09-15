clear;c;ca;
D = dir('*Mech.mat');
stdRange = [1:1000];
for ii = 1:length(D)
    tic
    ii
    load(D(ii).name,'out*')
    
    out_gD.raster = rmfield(out_gD,'raster');
    out_gJD.raster = rmfield(out_gJD,'raster');
    out_gND.raster = rmfield(out_gND,'raster');
    out_mD.raster = rmfield(out_mD,'raster');
    out_mJD.raster = rmfield(out_mJD,'raster');
    out_mND.raster = rmfield(out_mND,'raster');
    
    tic
    smSpike = zeros(length(stdRange),length(out_gD.tested_spikes));
    for jj = 1:length(stdRange)
        smSpike(jj,:) =smoothts(out_gD.tested_spikes','g',length(out_gD.tested_spikes),stdRange(jj));
    end
    smSpike = smSpike';
    
    timer = toc;
    fprintf('It took %.4f seconds to calculate the gaussian smoohthing \n%i samples and %i spikes\n',timer,length(out_gD.tested_spikes),sum(out_gD.tested_spikes))
    tic
    r(ii).predictiveG = corr(out_gD.predictiveModeY,smSpike);
    r(ii).predictiveGND = corr(out_gND.predictiveModeY,smSpike);
    r(ii).predictiveGJD = corr(out_gJD.predictiveModeY,smSpike);
    
    r(ii).predictiveM = corr(out_mD.predictiveModeY,smSpike);
    r(ii).predictiveMND = corr(out_mND.predictiveModeY,smSpike);
    r(ii).predictiveMJD = corr(out_mJD.predictiveModeY,smSpike);
    
    r(ii).generativeG = corr(out_gD.P,smSpike);
    r(ii).generativeGND = corr(out_gND.P,smSpike);
    r(ii).generativeGJD = corr(out_gJD.P,smSpike);
    
    r(ii).generativeM = corr(out_mD.P,smSpike);
    r(ii).generativeMND = corr(out_mND.P,smSpike);
    r(ii).generativeMJD = corr(out_mJD.P,smSpike);
    
    r(ii).noHistG = corr(out_gD.Y_noHist,smSpike);
    r(ii).noHistGND = corr(out_gND.Y_noHist,smSpike);
    r(ii).noHistGJD = corr(out_gJD.Y_noHist,smSpike);
    
    r(ii).noHistM = corr(out_mD.Y_noHist,smSpike);
    r(ii).noHistMND = corr(out_mND.Y_noHist,smSpike);
    r(ii).noHistMJD = corr(out_mJD.Y_noHist,smSpike);
    timer = toc;
    fprintf('It took %.4f seconds to perform correlations\n',timer)
        %%
        trace(ii).predictiveG = out_gD.predictiveModeY;
        trace(ii).predictiveGND = out_gND.predictiveModeY;
        trace(ii).predictiveGJD = out_gJD.predictiveModeY;
    
        trace(ii).predictiveM = out_mD.predictiveModeY;
        trace(ii).predictiveMND = out_mND.predictiveModeY;
        trace(ii).predictiveMJD = out_mJD.predictiveModeY;
    
        trace(ii).generativeG = out_gD.P;
        trace(ii).generativeGND = out_gND.P;
        trace(ii).generativeGJD = out_gJD.P;
    
        trace(ii).generativeM = out_mD.P;
        trace(ii).generativeMND = out_mND.P;
        trace(ii).generativeMJD = out_mJD.P;
    
        trace(ii).noHistG = out_gD.Y_noHist;
        trace(ii).noHistGND = out_gND.Y_noHist;
        trace(ii).noHistGJD = out_gJD.Y_noHist;
    
        trace(ii).noHistM = out_mD.Y_noHist;
        trace(ii).noHistMND = out_mND.Y_noHist;
        trace(ii).noHistMJD = out_mJD.Y_noHist;
        spikes{ii} = out_gD.tested_spikes;
    
    
end
save('Summary.mat','r','stdRange','trace','spikes')

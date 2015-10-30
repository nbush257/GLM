D = dir('*cleaned.mat');
for dd = 1:length(D)
    
    clearvars -except dd D keep
    if ~keep(dd)
        continue
    end
    dd
    D = dir('*cleaned.mat');
    load(D(dd).name)
    GLM_VG3D;
    fOut = [D(dd).name(1:end-12) '_RESULTS_noCV.mat'];
    clear D ydum
%     [x,y,~,AUC] = perfcurve(yTest,yOut,1);
    save(fOut)
end
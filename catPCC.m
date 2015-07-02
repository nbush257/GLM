dI = dir('*simGLM.mat')
allRG = NaN(length(dI),1);
allRM = NaN(length(dI),1);
allRG_med = NaN(length(dI),1);
allRG_dis = NaN(length(dI),1);
allRM_med = NaN(length(dI),1);
allRM_dis = NaN(length(dI),1);

for dIter = 1:length(dI)
    load(dI(dIter).name,'rG*','rM*')
    ca
    allRG(dIter) = rG;
    allRM(dIter) = rM;
    if exist('rG_med');
        allRG_med(dIter) = rG_med;
        allRG_dis(dIter) = rG_dis;
        allRM_med(dIter) = rM_med;
        allRM_dis(dIter) = rM_dis;
    end
    clearvars -except allRG allRM allRG_med allRG_dis allRM_med allRM_dis dIter dI
end
    
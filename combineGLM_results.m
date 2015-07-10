% get all vars
d = dir('*.mat')
medS = struct;
proxS = struct;
disS = struct;
for dd = 1:length(d)
    load(d(dd).name)
    geo(dd) = out_gD;
    mech(dd) = out_mD;
    g_No_Deriv(dd) = out_gND;
    m_No_Deriv(dd) = out_mND;
    g_JustDeriv(dd) = out_gJD;
    m_justDeriv(dd) = out_mJD;
    if exist('med','var')
    allMed{dd} = med;
    allProx{dd} = prox;
    allDis{dd} = dis;
    end
    
end

    



%compare 

f1 = figure 

bar([[out_gD.AUC out_gJD.AUC out_gND.AUC]' [out_mD.AUC out_mJD.AUC out_mND.AUC]'])
ax = gca;
ax.XTick = [1 2 3]
set(ax,'XTickLabel',{'Vars and Derivatives','Just derivatives','Just Variables'})
ax.Ylabel = 'AUC'
legend({'geometry','mechanics'});


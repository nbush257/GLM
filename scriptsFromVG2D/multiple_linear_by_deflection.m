for ii = 2:40
    tR = R{ii};
    tFR  =FR{ii};
    tM = M{ii};
    tFY = FY{ii};
    tFX = FX{ii};
    tTH = TH_o{ii};
    
    
    
    geo = stepwiselm([tR' tTH'],tFR','Upper','linear');
    useR(ii) = geo.VariableInfo.InModel(1);
    useTH(ii) = geo.VariableInfo.InModel(2);
        eval(['geo_' num2str(ii) ' = geo']);

    mech = stepwiselm([tFX' tM'],tFR','Upper','linear');
    useFX(ii) = mech.VariableInfo.InModel(1);
    useM(ii) = mech.VariableInfo.InModel(2);
            eval(['mech_' num2str(ii) ' = mech']);

    all = stepwiselm([tFX' tM' tR' tTH'],tFR','Upper','linear');
    eval(['all_' num2str(ii) ' = all']);
    useFX_all(ii) = all.VariableInfo.InModel(1);
    useM_all(ii) = all.VariableInfo.InModel(2);
    useR_all(ii) = all.VariableInfo.InModel(3);
    useTH_all(ii) = all.VariableInfo.InModel(4);
end

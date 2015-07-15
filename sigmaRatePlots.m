for ii = 1:length(names)
    if any(regexp(names(ii).name,'Rat1*'))
        nTitle{ii} = names(ii).name(9:12);
    else
        nTitle{ii} = names(ii).name([9:10 21 22 28:33])
    end
end

for ii = 1:length(r)
    subplot(6,7,ii)
    
    plot(stdRange(1:150),r(ii).noHistG(1:150));ln2
    axy(0,1)
    axx(0,150)
    title(nTitle{ii},'Interpreter','None')
    ho
    [gMax,gIdx] = max(r(ii).noHistG(1:150));
    if gIdx>=150
        gIdx = knee_pt(r(ii).noHistG(1:150));
        gMax = r(ii).noHistG(gIdx);
    end
   % plot(gIdx,gMax,'*')
    
    plot(stdRange(1:150),r(ii).noHistM(1:150));ln2
    [mMax,mIdx] = max(r(ii).noHistM(1:150));
    if mIdx>150
        mIdx = knee_pt(r(ii).noHistM(1:150));
        mMax = r(ii).noHistM(mIdx);
    end
    xlabel('\sigma (ms)')
    %ylabel('R^2')
    %plot(mIdx,mMax,'*')
end
for ii = 1:length(r)
    [~,sigmaStar(ii).generativeG] = max(r(ii).generativeG);
    [~,sigmaStar(ii).generativeGND] = max(r(ii).generativeGND);
    [~,sigmaStar(ii).generativeGJD] = max(r(ii).generativeGJD);
    [~,sigmaStar(ii).generativeM] = max(r(ii).generativeM);
    [~,sigmaStar(ii).generativeMND] = max(r(ii).generativeMND);
    [~,sigmaStar(ii).generativeMJD] = max(r(ii).generativeMJD);
    
    [~,sigmaStar(ii).predictiveG] = max(r(ii).predictiveG);
    [~,sigmaStar(ii).predictiveGND] = max(r(ii).predictiveGND);
    [~,sigmaStar(ii).predictiveGJD] = max(r(ii).predictiveGJD);
    [~,sigmaStar(ii).predictiveM] = max(r(ii).predictiveM);
    [~,sigmaStar(ii).predictiveMND] = max(r(ii).predictiveMND);
    [~,sigmaStar(ii).predictiveMJD] = max(r(ii).predictiveMJD);
    
    [~,sigmaStar(ii).noHistG] = max(r(ii).noHistG);
    [~,sigmaStar(ii).noHistGND] = max(r(ii).noHistGND);
    [~,sigmaStar(ii).noHistGJD] = max(r(ii).noHistGJD);
    [~,sigmaStar(ii).noHistM] = max(r(ii).noHistM);
    [~,sigmaStar(ii).noHistMND] = max(r(ii).noHistMND);
    [~,sigmaStar(ii).noHistMJD] = max(r(ii).noHistMJD);
    
    %%
    maxR(ii).generativeG =  max(r(ii).generativeG);
    maxR(ii).generativeGND =  max(r(ii).generativeGND);
    maxR(ii).generativeGJD =  max(r(ii).generativeGJD);
    maxR(ii).generativeM =  max(r(ii).generativeM);
    maxR(ii).generativeMND =  max(r(ii).generativeMND);
    maxR(ii).generativeMJD =  max(r(ii).generativeMJD);
    
    maxR(ii).predictiveG =  max(r(ii).predictiveG);
    maxR(ii).predictiveGND =  max(r(ii).predictiveGND);
    maxR(ii).predictiveGJD =  max(r(ii).predictiveGJD);
    maxR(ii).predictiveM =  max(r(ii).predictiveM);
    maxR(ii).predictiveMND =  max(r(ii).predictiveMND);
    maxR(ii).predictiveMJD =  max(r(ii).predictiveMJD);
    
    maxR(ii).noHistG =  max(r(ii).noHistG);
    maxR(ii).noHistGND =  max(r(ii).noHistGND);
    maxR(ii).noHistGJD =  max(r(ii).noHistGJD);
    maxR(ii).noHistM =  max(r(ii).noHistM);
    maxR(ii).noHistMND =  max(r(ii).noHistMND);
    maxR(ii).noHistMJD =  max(r(ii).noHistMJD);
end
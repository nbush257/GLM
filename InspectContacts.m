% d = dir('*cell_1khz.mat');
% 
% for ii = 4
%     % Be paranoid -- clear the data
%     clearvars -except d ii
%     load(d(ii).name);
%     d(ii).name
    
    mechdata  =  mech_85_justRCCR(:,2:end)';
    geodata  =   geo_85_justRCCR(:,2:end)';
    spikevec = mech_85_justRCCR(:,1)';
    st = find(mech_85_justRCCR(:,1)==1)';
    
    FX_noFilt = mech_noFilt_justRCCR(:,2)';
    FX_25 = mech_25_justRCCR(:,2)';
    FX_85 = mech_85_justRCCR(:,2)';
    
    FY_noFilt = mech_noFilt_justRCCR(:,3)';
    FY_25 = mech_25_justRCCR(:,3)';
    FY_85 = mech_85_justRCCR(:,3)';
    
    M_noFilt = mech_noFilt_justRCCR(:,4)';
    M_25 = mech_25_justRCCR(:,4)';
    M_85 = mech_85_justRCCR(:,4)';
    
    R_noFilt = geo_noFilt_justRCCR(:,2)';
    R_25 = geo_25_justRCCR(:,2)';
    R_85 = geo_85_justRCCR(:,2)';
    
    TH_noFilt = geo_noFilt_justRCCR(:,3)';
    TH_25 = geo_25_justRCCR(:,3)';
    TH_85 = geo_85_justRCCR(:,3)';
    
    
    
    
    clf;
    [C,NC,c,nc,starts,stops,contactDurations,ubmech,ubgeo] = FindContacts(mechdata, geodata, 50, st);
    for jj = 1:length(starts)
        text(mean([starts(jj) stops(jj)]),1,num2str(jj),'FontSize',16,'HorizontalAlignment','center')
    end
    rm  = input('Which contacts to remove?');
    modify = input('Which contacts to modify?');
    for jj = 1:length(rm)
        FX_noFilt(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        FX_25(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        FX_85(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        
        FY_noFilt(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        FY_25(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        FY_85(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        
        M_noFilt(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        M_25(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        M_85(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        
        R_noFilt(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        R_25(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        R_85(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        
        TH_noFilt(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        TH_25(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        TH_85(:,starts(rm(jj)):stops(rm(jj))) = NaN;
        
        C(starts(rm(jj)):stops(rm(jj))) = 0;
        
    end
    
    
    toNaN = [];
    startToNaN = [];
    stopToNaN = []
    for jj = 1:length(modify)
        axx(starts(modify(jj))-10,stops(modify(jj))+10);
        toNaN = ginput;
        toNaN = toNaN(:,1);
        startToNaN = [startToNaN;toNaN(1:2:end)];
        stopToNaN = [stopToNaN; toNaN(2:2:end)];
    end
    startToNaN = round(startToNaN);
    stopToNaN = round(stopToNaN);
    
    
    for jj = 1:length(startToNaN)
        mechdata(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        geodata(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        C(startToNaN(jj):stopToNaN(jj)) = NaN;
        
        FX_noFilt(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        FX_25(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        FX_85(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        
        FY_noFilt(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        FY_25(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        FY_85(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        
        M_noFilt(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        M_25(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        M_85(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        
        R_noFilt(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        R_25(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        R_85(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        
        TH_noFilt(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        TH_25(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        TH_85(:,startToNaN(jj):stopToNaN(jj)) = NaN;
        
        C(startToNaN(jj):stopToNaN(jj)) = 0;
        
    end
    
    starts = find(diff(C)==1)+1;
    stops = find(diff(C)==-1);
    for jj = 1:length(starts)
        c{jj} = starts(jj):stops(jj);
    end
    
%     
%     save([d(ii).name(1:end-4) '_cleaned.mat'],'FX_noFilt','FX_25','FX_85','FY_noFilt','FY_25','FY_85','M_noFilt','M_25','M_85','R_25','R_noFilt','R_85','TH_25','TH_85','TH_noFilt','C','c','starts','stops','st','spikevec')
%     
%     
% end;

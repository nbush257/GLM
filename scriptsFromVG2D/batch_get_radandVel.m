dire = dir('rat2015_15*')
ca
for ddire = 1:length(dire)
        clearvars -except dire ddire

    load(dire(ddire).name)
    
    
    Geo.vel = Geo.filtAll(:,2);
    Geo.vel = cdiff(Geo.vel);
    
    Geo.vel = interpNaNFilt(Geo.vel,1000,85);
    Geo.vel(isnan(Geo.vel)) = 0;
    ca
    save(dire(ddire).name,'Geo','Mech','spikevec','dis','prox','med','geo_*','mech_*')
end
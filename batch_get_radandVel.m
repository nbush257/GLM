d = dir('rat20*')
ca
for dd = 12:length(d)
        clearvars -except d dd

    load(d(dd).name)
    [prox,med,dis] = getRadialDistanceGroup(Geo);
    prox = logical(prox);
    med = logical(med);
    
    Geo.vel = Geo.filtTH;
    Geo.vel = cdiff(Geo.vel);
    
    Geo.vel = interpNaNFilt(Geo.vel,1000,85);
    Geo.vel(isnan(Geo.vel)) = 0;
    ca
    save(d(dd).name,'Geo','Mech','spikevec','dis','prox','med')
end
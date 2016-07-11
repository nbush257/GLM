d = dir('*concat2.mat');

for ii = 1:length(d)
    load(d(ii).name)
    FX = medfilt1(FX);
    FY = medfilt1(FY);
    M = medfilt1(M);
    TH = medfilt1(TH);
    TH(TH==0) = NaN;
    R = medfilt1(R);
    V = cdiff(TH);
    V = interpNaNFilt(V,1000,85);
    C(isnan(C)) = 0;
    
    
    Mech.filtAll = [FX FY M];
    Mech.filtFX = FX;
    Mech.filtFY = FY;
    Mech.filtM = M;
    
    Geo.filtAll = [R TH V];
    Geo.filtR = R;
    Geo.filtTH = TH;
    Geo.filtV = V;
    
    C = logical(C);
    
    save([d(ii).name(1:end-11) 'toGLM.mat'],'Geo','Mech','spikevec')
    
    
    
end
    

for ii = 1:length(varis)
        subplot(211)
        scatter(varis(ii).TH,varis(ii).FX,10*ones(size(varis(ii).FX)),varis(ii).R)
        colorbar
        xlabel('Angular displcement(deg)')
        ylabel('FX (N)')
        title(num2str(ii))
        subplot(212)
        scatter(varis(ii).TH,varis(ii).M,10*ones(size(varis(ii).FX)),varis(ii).R)
        colorbar
        xlabel('Angular displcement(deg)')
        ylabel('M (N-m)')
    pause
end


        
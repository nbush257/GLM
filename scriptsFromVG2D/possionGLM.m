% poisson
clear
load('rat2015_15_JUN11_VG_B2_t01_cell_1_toGLM.mat')
C = logical(C);
C(1) = 0;
C(end) = 0;

if size(mech_85,2)>size(mech_85,1)
    mech_85 = mech_85';
    geo_85 = geo_85';
end

if sum(spikevec)<50
    continue
end
if isrow(spikevec);spikevec = spikevec'; end;
if isrow(C);C = C'; end;

if noProx & exist('prox','var')
    C(prox) = [];
    mech_85(prox,:) = [];
    geo_85(prox,:) = [];
    spikevec(prox) = [];
    med(prox) = [];
    dis(prox) = [];
end




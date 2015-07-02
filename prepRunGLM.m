%%
geo = geo_85;
mech = mech_85;
%%
C(isnan(C)) = 0;
C = logical(C);
blipsUp = strfind(C,[0 1 0])+1;
blipsDown = strfind(C,[1 0 1])+1;
C(blipsUp) = 0;
C(blipsDown) = 1;
spikevec(spikevec>1)=1;

%
%% fix TH
TH = geo(2,:);

start = find(diff(C)==1)+1;
stop = find(diff(C)==-1);


THnew = nan(size(TH));
for ii = 1:length(start)
    non_nan = find(~isnan(TH(start(ii):stop(ii))));
    if isempty(non_nan)
        continue
    end
    THnew(start(ii):stop(ii)) = TH(start(ii):stop(ii))-TH(start(ii)+non_nan(1)-1);
end
geo(2,:) = THnew;
geo(1,geo(1,:)==0) = NaN;
%% prep geometric
geo(:,1) = 0;
geo(:,end) = 0;
geo(2,:) = naninterp(geo(2,:));
geo(:,~C)=0;
start = find(diff(C)==1)+1;
stop = find(diff(C)==-1);
newSpike = spikevec;
newGeo = [];newSpike = [];
if C(1) ==1
    start = [1 start];
end
if C(end)==1
    stop = [stop length(C)];
end
if iscolumn(spikevec);spikevec = spikevec';end
for ii = 1:length(start)
    temp = [[nan(1,15);zeros(1,15)] geo(:,start(ii):stop(ii)) [nan(1,15);zeros(1,15)]];
    tempSpike = [zeros(1,15) spikevec(start(ii):stop(ii)) zeros(1,15)];
    newGeo  = [newGeo temp];
    newSpike = [newSpike tempSpike];
end
if isempty(start)
    newGeo = geo;
end
%% Mechanical

mech(:,1) = 0;
mech(:,end) = 0;
mech = naninterp(mech);
mech(:,~C)=0;

newMech = [];newSpike = [];
% 
% for ii = 1:length(start)
%     temp = [[zeros(3,15)] mech(:,start(ii):stop(ii)) [zeros(3,15)]];
%     tempSpike = [zeros(1,15) spikevec(start(ii):stop(ii)) zeros(1,15)];
%     newMech  = [newMech temp];
%     newSpike = [newSpike tempSpike];
% end
if isempty(start)
    newGeo = geo;
    newMech = mech;
    newSpike = spikevec;
end
geo = newGeo;
mech = newMech;
spikes = newSpike;

%% Radial distance metric
toKeep = [];
for ii = 1:length(start)
    toKeep = [toKeep start(ii)-15:stop(ii)+15];
end
dis = dis(toKeep);
prox = prox(toKeep);
med = med(toKeep);


function [prox,med,dis] = getRadialDistanceGroup(Geo)
%% function [prox,med,dis] = getRadialDistanceGroup(Geo)
% Takes the Radial distance from the Geo struct and builds a GUI to
% determine which sections are proximal, medial, and distal. Returns
% logical vectors te length of Geo.R


%% input handling
if ~isstruct(Geo)
    error('input is not a structure')
end

h = figure;
%% Find Proximal
plot(Geo.filtR)
title('Find Most proximal, right click to exit')
but = 0;
x = [];
count = 0;

prox = zeros(size(Geo.filtR));

while but~=3 % right click exits
    count = count+1;
    zoom on
    pause
    
    X = xlim;
    Y = ylim;
    [x,~,but] = ginput(2);
    
    if isempty(x)
        but = 3;
    continue
    end
    if x(2)>length(Geo.filtR)
        x(2) = length(Geo.filtR);
    end
    
    if x(1)<1
        x(1) = 1;
    end
    prox(x(1):x(2)) = 1;
    cla
    plot(Geo.filtR)
    ho

    plot(find(prox),Geo.filtR(find(prox)),'g')
    
    axx(X(1),X(2));
    axy(Y(1),Y(2));
end

%% Find Distal
dis = zeros(size(Geo.filtR));
title('Find Most distal, right click to exit')
but = 0;
while but~=3 % right click to exit
    count = count+1;
    zoom on
    pause
    X = xlim;
    Y = ylim;
    [x,~,but] = ginput(2);
    if isempty(x)
        but = 3;
    continue
    end
    
    if x(1)<1
        x(1) = 1;
    end
    
    
    if x(2)>length(Geo.filtR)
        x(2) = length(Geo.filtR);
    end
    dis(x(1):x(2)) = 1;
    cla
    plot(Geo.filtR)
    ho
    plot(find(prox),Geo.filtR(find(prox)),'g')
    plot(find(dis),Geo.filtR(find(dis)),'r')
    axx(X(1),X(2));
    axy(Y(1),Y(2));
    
    
end

%% get medial as exclusion of other distances
prox = logical(prox);
dis = logical(dis);
med = ~isnan(Geo.R) & ~prox & ~dis;
plot(Geo.R)
title('Find Most proximal, right click to exit')
but = 0;
x = [];
count = 0;

prox = zeros(size(Geo.R));

while but~=3
    count = count+1;
    zoom on
    pause
    [x,~,but] = ginput(2);
    prox(x(1):x(2)) = 1;
    ho
    plot(find(prox),Geo.R(find(prox)),'g')
end


dis = zeros(size(Geo.R));

title('Find Most distal, right click to exit')
but = 0;
while but~=3
    count = count+1;
    zoom on
    pause
    [x,~,but] = ginput(2);
    dis(x(1):x(2)) = 1;
    ho
    plot(find(dis),Geo.R(find(dis)),'r')
end


med = C & ~prox & ~dis;
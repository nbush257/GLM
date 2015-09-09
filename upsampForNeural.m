function xOut = upsampForNeural(x,y,frameCapSamps,startFrame,endFrame,numFrames)
% this allows you to upsample a vector x to fit another vector y (used for
% comparing mechanics to neural)

x_upsamp = nan(size(y));

for ii = 1:numFrames-1
    if any(isnan(ceil(frameCapSamps(startFrame+ii-1)):ceil(frameCapSamps(startFrame+ii))))
        continue
    end
    x_upsamp(ceil(frameCapSamps(startFrame+ii-1)):ceil(frameCapSamps(startFrame+ii))) = x(ii);
end

x_pts = nan(size(y));

for ii = 1:numFrames-1 
    if isnan(ceil(frameCapSamps(startFrame+ii-1)))
        continue
    end
    x_pts(ceil(frameCapSamps(startFrame+ii-1))) = x(ii);
end

xOut=x_pts;

bd_x=isnan(x_pts) & ~isnan(x_upsamp);

gd_x=find(~bd_x);

bd_x([1:(min(gd_x)-1) (max(gd_x)+1):end])=0;

xOut(bd_x)=interp1(gd_x,x_pts(gd_x),find(bd_x));

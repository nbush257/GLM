function x_down = downsampForNeural(x,framesamps)
%% function x_down = downsampForNeural(x,framesamps)
% this function will downsample a mechanical or kinematic variable to the
% original framerate of the video. There are small differences between the
% original vector and the upsampled-downsampled vector, so do not convert
% frequently if it can be avoided. 
%%
x_down = nan(length(framesamps),size(x,2));
x_down = x(framesamps,:);

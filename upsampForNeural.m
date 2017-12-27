function x_upsamp = upsampForNeural(x,y,framesamps,varargin)
%% function x_upsamp = upsampForNeural(x,y,framesamps,[win_size])
% this allows you to upsample a vector x to fit another vector y (used for
% comparing mechanics to neural)
% INPUTS:       x - Mechanics matrix 
%               y - boolean vector of spikes sampled at acquisition
%                   sampling rate
%               framesamps - indices of trigger onsets for each frame in
%                  acquisition samples
%               [win_size] - a window over which to interpolate nans. Should
%               be at least greater than the inter-frame interval (in
%               samples) Default = 2*interframe interval
% 
% OUTPUTS:      x_upsamp - upsampled mechanics matrix
% NEB 2017_02_14
%% get interpolation window
if nargin == 3
    inter_frame_interval = max(diff(framesamps));
    assert(inter_frame_interval<1000,'The interframe interval is too large, something must be wrong')
    win_size = ceil(inter_frame_interval*2);
elseif nargin == 4
    win_size = varargin{1};
else
    error('Too many arguments')
end

%% init 

% check to see if matrix has time along length of column (i.e., each row is
% a timepoint)
if size(x,2)>size(x,1)
    error('Stimulus matrix is likely not oriented with time along length of column.')
end

    
x_upsamp = nan(size(y,1),size(x,2));

% removes the first frame sample indicator in the case that the first sampl
% looks like a frame
% if framesamps(1) == 1
%     framesamps(1) = [];
% end
assert(length(x) == length(framesamps),'Number of captured frames is not the same as the number of stimulus samples.')


%% align
x_upsamp(framesamps,:) = x;
% interpolate
for ii = 1:size(x_upsamp,2)
    x_upsamp(:,ii) = InterpolateOverNans(x_upsamp(:,ii),win_size);
end

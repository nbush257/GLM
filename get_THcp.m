function TH_cp = get_THcp(xs,ys,CP,C)
%% function TH_cp = get_THcp(xs,ys,CP,C)
% INPUTS:   xs - a N cell array of whisker x points (each cell is a frame)
%           ys - a N element cell array of whisker y points 
%           CP - a Nx2 matrix of contact points 
%           C -  a N x 1 logical indicating contact
% OUTPUTS: TH_cp - an N x 1 vector of the angle theta_cp: the angle between
% the whisker base tangent and the line connecting the basepoint with the
% contactpoint.
%==================================================
% takes x and y pts of a whisker, along with the contact point and a
% contact binary to calculate theta_cp which is the angle between:
% 
% 1) the tangent at the whisker base (found using a linear fit to the first 10% of
% the whisker) 
% 2) the line connecting the basepoint and the contact point.
%
% This angle is in whisker centered coordinates proper. This code only
% computes the angle during times of contact. Otherwise, it is nan.
% =================================================================
% Nick Bush 2016_05_11
%% handle inputs
% makes CP an n x 2 column matrix. If CP is a 2x2 matrix, don't flip
if size(CP,1)==2 && size(CP,1)~=size(CP,2)
    CP = CP';
end

% make C a column vector
if isrow(C); C = C';end

C = logical(C);

%% main loop
% init theta_cp
TH_cp = nan(size(C));
disp('Calculating TH_cp...')
warning('off')
for ii = 1:length(xs) %iterate over every frame
    
    % skip if no contact
    if ~C(ii)
        continue
    end
    % skip if there is no whisker
    if isempty(xs{ii})
        continue
    end
    % skip if the whisker is too small
    l = length(xs{ii});
    if round(l/10)==0
        continue
    end
    
    % get basepoint
    x1 = xs{ii}(1);
    y1 = ys{ii}(1);
    
    % linear fit to first 10% of whisker to get the tangent.
    p = polyfit(xs{ii}(1:round(l/10)),ys{ii}(1:round(l/10)),1);
    xq = xs{ii}(round(l/10));
    yq = polyval(p,xq);
    
    % Calculate angle of tangent wrt video
    TH_linear = atan2(yq-y1,xq-x1)*180/pi;
    
    % Calc angle of ContactPoint wrt video
    TH = atan2(CP(ii,2)-y1,CP(ii,1)-x1)*180/pi;
    
    % subtract base tangent from contactpoint
    TH_cp(ii) = TH-TH_linear;
    
end
warning('on')
disp('wrapping')
% Wrap TH_cp
TH_cp(TH_cp>nanmedian(TH_cp)+180) = TH_cp(TH_cp>nanmedian(TH_cp)+180)-360;
TH_cp(TH_cp<nanmedian(TH_cp)-180) = TH_cp(TH_cp<nanmedian(TH_cp)-180)+360;
disp('done')

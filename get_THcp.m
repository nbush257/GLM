function TH_cp = get_THcp(xs,ys,CP,C)
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

% Wrap TH_cp
TH_cp(TH_cp>nanmedian(TH_cp)+180) = TH_cp(TH_cp>nanmedian(TH_cp)+180)-360;
TH_cp(TH_cp<nanmedian(TH_cp)-180) = TH_cp(TH_cp<nanmedian(TH_cp)-180)+360;

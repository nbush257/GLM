function [X_upsamp,spbool,spt_upsamp,tvec] = resamp(X,spt,sr,frametimes,new_sr)
%% function [X,spbool,tvec] = resamp(X,spt,sr,frametimes,new_sr)
% outputs a stimulus matrix as well as a spbool with each time bin the
% inverso of new sr (ms bins means new_sr = 1000)

maxt = frametimes(end); % in seconds
tvec = [0:maxt*new_sr];

frametimes_upsam = round(frametimes*new_sr);
X_upsamp = nan(length(tvec),size(X,2));
X_upsamp(frametimes_upsam,:) = X;
for ii = 1:size(X,2)
X_upsamp(:,ii) = InterpolateOverNans(X_upsamp(:,ii),sr/new_sr);
X_upsamp(:,ii) = naninterp(X_upsamp(:,ii));
end

spbool = zeros(size(tvec));
spt_upsamp = ceil(spt*new_sr);
u = unique(spt_upsamp);
n=histc(spt_upsamp,u);
spbool(u) = n;



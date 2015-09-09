%%
frameCapTimes(frameCapTimes>time(length(time)))=NaN;
timeSubSamp = floor(time*1000);
timeMs = 1:timeSubSamp(end);
frameCapTimes_ms = frameCapTimes*1000;


numFrames = endFrame - startFrame +1;

f = bwfilt(analogData.channel_0,sr,300,6000);
RCCR_start_samp = frameCapSamps(RCCR(1)+startFrame);
RCCR_end_samp = frameCapSamps(RCCR(2)+startFrame);

RCCR_start_ms = timeSubSamp(RCCR_start_samp);
RCCR_end_ms = timeSubSamp(RCCR_end_samp);

startSamp = frameCapSamps(startFrame);
endSamp = frameCapSamps(endFrame);
%%
neural_logical = times(cell).times>startSamp & times(cell).times<endSamp;% pick your cell
neural_clip = times(cell).times(neural_logical); % pick your cell
neural_ms = timeSubSamp(neural_clip);
[~,idx] = unique(neural_ms);
double = setdiff(1:length(neural_ms),idx);
neuralWord = zeros(size(timeMs));
neuralWord(unique(neural_ms)) = 1;
for ii = 1:length(double)
    neuralWord(neural_ms(double(ii))) = neuralWord(neural_ms(double(ii)))+1;
end


plot(f);
ho
plot(times(cell).times,f(times(cell).times),'o')
plot(neural_clip,f(neural_clip),'o')

title('make sure this is the cell you want')
pause
close all force
pause(.01)
%% Remove large outliers

FX(abs(FX)>.1) = NaN;
FY(abs(FY)>.1) = NaN;
M(abs(M)>.1) = NaN;



plot(FX)
zoom on;
title('zoom')
pause;
title('click above, then below')
thresh = ginput(2);
if numel(thresh)==4
FX(FX>thresh(1,2))=NaN;
FX(FX<thresh(2,2))=NaN;
end


plot(FY)
zoom on;
title('zoom')
pause;
title('click above, then below')
thresh = ginput(2);

if numel(thresh)==4
FY(FY>thresh(1,2))=NaN;
FY(FY<thresh(2,2))=NaN;
end
plot(M)
zoom on;
title('zoom')
pause;
title('click above, then below')
thresh = ginput(2);
if numel(thresh)==4
M(M>thresh(1,2))=NaN;
M(M<thresh(2,2))=NaN;
end
ca
pause(.01)
%% Get R
for ii = 1:numFrames
    if isempty(xs{ii}) | isnan(CP(ii,1))
        R(ii) = NaN;
        continue
    end
    R(ii) = sqrt((xs{ii}(1)-CP(ii,1)).^2 + (ys{ii}(1)-CP(ii,2)).^2);
end

%% Wrap TH

% getTH
TH_CP = nan(size(xs));
TH = nan(size(xs));

for ii = 1:length(xs)
    if isempty(xs{ii}) 
        continue
    end
    x1 = xs{ii}(1);
    y1 = ys{ii}(1);
    l = length(xs{ii});
    if round(l/10)==0
        continue
    end
    ye = ys{ii}(round(l/10));
    xe = xs{ii}(round(l/10));
    
    TH_linear(ii) = atan2(ye-y1,xe-x1)*180/pi;
    if C(ii)
        TH_CP(ii) = atan2(CP(ii,2)-y1,CP(ii,1)-x1)*180/pi;
        TH(ii) = TH_CP(ii)-TH_linear(ii);
    end
end

TH(TH>nanmean(TH)+180) = TH(TH>nanmean(TH)+180)-360;
TH(TH<nanmean(TH)-180) = TH(TH<nanmean(TH)-180)+360;

%% Lucie NaN Filter

% FX

FX_85 = interpNaNFilt(FX,300,85);
FX_25 = interpNaNFilt(FX,300,25);
% FY
FY_85 = interpNaNFilt(FY,300,85);
FY_25 = interpNaNFilt(FY,300,25);

% M
M_85 = interpNaNFilt(M,300,85);
M_25 = interpNaNFilt(M,300,25);

%R 
if isrow(R);R = R';end
R_85 = interpNaNFilt(R,300,85);
R_25 = interpNaNFilt(R,300,25);

% TH
if isrow(TH);TH = TH'; end;
TH_85 = interpNaNFilt(TH,300,85);
TH_25 = interpNaNFilt(TH,300,25);

%% Interpolate
FX_upsamp = upsampForNeural(FX,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);
FX_85_upsamp = upsampForNeural(FX_85,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);
FX_25_upsamp = upsampForNeural(FX_25,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);

FY_upsamp = upsampForNeural(FY,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);
FY_85_upsamp = upsampForNeural(FY_85,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);
FY_25_upsamp = upsampForNeural(FY_25,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);

M_upsamp = upsampForNeural(M,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);
M_85_upsamp = upsampForNeural(M_85,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);
M_25_upsamp = upsampForNeural(M_25,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);

R_upsamp = upsampForNeural(R,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);
R_85_upsamp = upsampForNeural(R_85,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);
R_25_upsamp = upsampForNeural(R_25,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);

TH_upsamp = upsampForNeural(TH,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);
TH_85_upsamp = upsampForNeural(TH_85,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);
TH_25_upsamp = upsampForNeural(TH_25,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);

C_upsamp = upsampForNeural(C,timeMs,frameCapTimes_ms,startFrame,endFrame,numFrames);


%% Save all data

mechMat_noFilt = [neuralWord' FX_upsamp' FY_upsamp' M_upsamp'];
mechMat_85 = [neuralWord' FX_85_upsamp' FY_85_upsamp' M_85_upsamp'];
mechMat_25 = [neuralWord' FX_25_upsamp' FY_25_upsamp' M_25_upsamp'];
geoMat_noFilt = [neuralWord' R_upsamp' TH_upsamp'];
geoMat_85 = [neuralWord' R_85_upsamp' TH_85_upsamp'];
geoMat_25 = [neuralWord' R_25_upsamp' TH_25_upsamp'];



%save just this clip
clipMechMat_noFilt = nan(size(mechMat_noFilt));
clipMechMat_85 = nan(size(mechMat_85));
clipMechMat_25 = nan(size(mechMat_25));
clipGeoMat_noFilt = nan(size(geoMat_noFilt));
clipGeoMat_85 = nan(size(geoMat_85));
clipGeoMat_25 = nan(size(geoMat_25));


clipMechMat_noFilt(RCCR_start_ms:RCCR_end_ms,:) = mechMat_noFilt(RCCR_start_ms:RCCR_end_ms,:);
clipMechMat_85(RCCR_start_ms:RCCR_end_ms,:) = mechMat_85(RCCR_start_ms:RCCR_end_ms,:);
clipMechMat_25(RCCR_start_ms:RCCR_end_ms,:) = mechMat_25(RCCR_start_ms:RCCR_end_ms,:);
clipGeoMat_noFilt(RCCR_start_ms:RCCR_end_ms,:) = geoMat_noFilt(RCCR_start_ms:RCCR_end_ms,:);
clipGeoMat_85(RCCR_start_ms:RCCR_end_ms,:) = geoMat_85(RCCR_start_ms:RCCR_end_ms,:);
clipGeoMat_25(RCCR_start_ms:RCCR_end_ms,:) = geoMat_25(RCCR_start_ms:RCCR_end_ms,:);
C_upsamp = C_upsamp(RCCR_start_ms:RCCR_end_ms);
%% Get justRCCR

notRCCR = all(isnan(clipMechMat_noFilt),2);

mech_noFilt_justRCCR = clipMechMat_noFilt(~notRCCR,:);
mech_85_justRCCR = clipMechMat_85(~notRCCR,:);
mech_25_justRCCR = clipMechMat_25(~notRCCR,:);

geo_noFilt_justRCCR = clipGeoMat_noFilt(~notRCCR,:);
geo_85_justRCCR = clipGeoMat_85(~notRCCR,:);
geo_25_justRCCR = clipGeoMat_25(~notRCCR,:);
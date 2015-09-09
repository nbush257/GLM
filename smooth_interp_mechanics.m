% TAG_Mech = 'rat2015_06_FEB26_VG_C1_t01_F060001F069590';
% TAG_E2D = 'Vg_DATA_PROC_rat2015_06_FEB26_vg_C1_t01_Top_F060001F069590_E2D';
% load(['C:\Users\nbush257\Documents\hartmann_lab\data\for_E2D_emergency\forJen\' TAG_Mech '_2DMech_w_neural.mat'])
% load(['C:\Users\nbush257\Documents\hartmann_lab\data\for_E2D_emergency\E2D_results\' TAG_E2D ])


%% Get R
for ii = 1:numFrames
    if isempty(xs{ii}) | isnan(CP(ii,1))
        R(ii) = NaN;
        continue
    end
    R(ii) = sqrt((xs{ii}(1)-CP(ii,1)).^2 + (ys{ii}(1)-CP(ii,2)).^2);
end
%%
TH(TH>nanmean(TH)+180) = TH(TH>nanmean(TH)+180)-360;
TH(TH<nanmean(TH)-180) = TH(TH<nanmean(TH)-180)+360;

%%
FX(abs(FX)>.001) = NaN;
FY(abs(FY)>.001) = NaN;
M(abs(M)>.002) = NaN;

FX_raw = FX;
FY_raw = FY;
M_raw = M;

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
FX_upsamp = upsampForNeural(FX,f,frameCapSamps,startFrame,endFrame,numFrames);
FX_85_upsamp = upsampForNeural(FX_85,f,frameCapSamps,startFrame,endFrame,numFrames);
FX_25_upsamp = upsampForNeural(FX_25,f,frameCapSamps,startFrame,endFrame,numFrames);

FY_upsamp = upsampForNeural(FY,f,frameCapSamps,startFrame,endFrame,numFrames);
FY_85_upsamp = upsampForNeural(FY_85,f,frameCapSamps,startFrame,endFrame,numFrames);
FY_25_upsamp = upsampForNeural(FY_25,f,frameCapSamps,startFrame,endFrame,numFrames);

M_upsamp = upsampForNeural(M,f,frameCapSamps,startFrame,endFrame,numFrames);
M_85_upsamp = upsampForNeural(M_85,f,frameCapSamps,startFrame,endFrame,numFrames);
M_25_upsamp = upsampForNeural(M_25,f,frameCapSamps,startFrame,endFrame,numFrames);

R_upsamp = upsampForNeural(R,f,frameCapSamps,startFrame,endFrame,numFrames);
R_85_upsamp = upsampForNeural(R_85,f,frameCapSamps,startFrame,endFrame,numFrames);
R_25_upsamp = upsampForNeural(R_25,f,frameCapSamps,startFrame,endFrame,numFrames);

TH_upsamp = upsampForNeural(TH,f,frameCapSamps,startFrame,endFrame,numFrames);
TH_85_upsamp = upsampForNeural(TH_85,f,frameCapSamps,startFrame,endFrame,numFrames);
TH_25_upsamp = upsampForNeural(TH_25,f,frameCapSamps,startFrame,endFrame,numFrames);

% Save all data
mechMat_noFilt = [neuralWord FX_upsamp FY_upsamp M_upsamp];
mechMat_85 = [neuralWord FX_85_upsamp FY_85_upsamp M_85_upsamp];
mechMat_25 = [neuralWord FX_25_upsamp FY_25_upsamp M_25_upsamp];
geoMat_noFilt = [neuralWord R_upsamp TH_upsamp];
geoMat_85 = [neuralWord R_85_upsamp TH_85_upsamp];
geoMat_25 = [neuralWord R_25_upsamp TH_25_upsamp];


% save just this clip
clipMechMat_noFilt = nan(size(mechMat_noFilt));
clipMechMat_85 = nan(size(mechMat_85));
clipMechMat_25 = nan(size(mechMat_25));
clipGeoMat_noFilt = nan(size(geoMat_noFilt));
clipGeoMat_85 = nan(size(geoMat_85));
clipGeoMat_25 = nan(size(geoMat_25));

clipMechMat_noFilt(RCCR_start_samp:RCCR_end_samp,:) = mechMat_noFilt(RCCR_start_samp:RCCR_end_samp,:);
clipMechMat_85(RCCR_start_samp:RCCR_end_samp,:) = mechMat_85(RCCR_start_samp:RCCR_end_samp,:);
clipMechMat_25(RCCR_start_samp:RCCR_end_samp,:) = mechMat_25(RCCR_start_samp:RCCR_end_samp,:);
clipGeoMat_noFilt(RCCR_start_samp:RCCR_end_samp,:) = geoMat_noFilt(RCCR_start_samp:RCCR_end_samp,:);
clipGeoMat_85(RCCR_start_samp:RCCR_end_samp,:) = geoMat_85(RCCR_start_samp:RCCR_end_samp,:);
clipGeoMat_25(RCCR_start_samp:RCCR_end_samp,:) = geoMat_25(RCCR_start_samp:RCCR_end_samp,:);
%% Get justRCCR

notRCCR = all(isnan(clipMechMat_noFilt),2);

mech_noFilt_justRCCR = clipMechMat_noFilt(~notRCCR,:);
mech_85_justRCCR = clipMechMat_85(~notRCCR,:);
mech_25_justRCCR = clipMechMat_25(~notRCCR,:);

geo_noFilt_justRCCR = clipGeoMat_noFilt(~notRCCR,:);
geo_85_justRCCR = clipGeoMat_85(~notRCCR,:);
geo_25_justRCCR = clipGeoMat_25(~notRCCR,:);

%%


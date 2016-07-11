%% alginNeuralAnd2D
% % Lots of preprocessing, then upsampling inputs.
% clear;ca
% [fName,pName] = uigetfile('*varConcat.mat','Load in the concatenated variables from E2D');
%
%
% oldir = pwd;
% cd P:\toGLM_V4\thetaFix\neural
% fName
% [fName2,pName2] = uigetfile('*.mat',['Load in the sorted neural data for ' fName]);
% fName2
% load([pName2 fName2]);
% load([pName fName]);
% cd(oldir);
% fOut = regexp(fName2,'_t\d{2}');fOut = fName2(1:fOut+3);

function AlignNeuralAndE2D(fName,fName2,useCells)
plotTGL = 0;
cInterpWin = 20;
times = [];
load(fName)
load(fName2)
fOut = regexp(fName2,'_t\d{2}');fOut = fName2(1:fOut+3);
%% Neural data to ms
frameCapTimes(frameCapTimes>time(length(time)))=NaN;
timeSubSamp = floor(time*1000);
timeMs = 1:timeSubSamp(end);
frameCapTimes_ms = frameCapTimes*1000;
if plotTGL
    figure
    for ii = 1:length(shapes)
        subplot(1,length(shapes),ii)
        plot(shapes(ii).shapes','k')
        title(num2str(ii))
        axy(-150,150)
    end
end
% useCells = input('Which cells do you want to use?');
if length(useCells) ==0
    useCells = 1:length(times);
end

neural_clip = [];
neural_ms = [];
neural_word = zeros(length(timeMs),length(useCells));
for ii = 1:length(useCells)
    neural_clip = times(useCells(ii)).times; % pick your cell
    neural_ms = timeSubSamp(neural_clip);
    neural_word(neural_ms,ii) = 1;
end


%% make vars easier to access
xs = E2D.xs;
ys = E2D.ys;
C = E2D.C;
CP = E2D.CP;
TH_CP = nan(size(xs));
TH = nan(size(xs));
%% Clean C variable
tempC = C;
tempC(tempC==0) = nan;
C = InterpolateOverNans(tempC,cInterpWin);
C(isnan(C)) = 0;
E2D.C = logical(C);
%% getTH_cp
for ii = 1:length(xs)
    if C(ii)
        if isempty(xs{ii})
            continue
        end
        x1 = xs{ii}(1);
        y1 = ys{ii}(1);
        
        l = length(xs{ii});
        if round(l/10)==0
            continue
        end
        
        %possible linear fit to this section of the whisker.
        ye = ys{ii}(round(l/10));
        xe = xs{ii}(round(l/10));
        
        TH_linear(ii) = atan2(ye-y1,xe-x1)*180/pi;
        
        TH_CP(ii) = atan2(CP(ii,2)-y1,CP(ii,1)-x1)*180/pi;
        
        TH(ii) = TH_CP(ii)-TH_linear(ii);
    end
end

% wrap
TH(TH>nanmedian(TH)+180) = TH(TH>nanmedian(TH)+180)-360;
TH(TH<nanmedian(TH)-180) = TH(TH<nanmedian(TH)-180)+360;
% output
if isrow(TH);TH = TH';end
E2D.TH_cp = TH;
E2D_medfilt.TH_cp = medfilt1(TH);
%% get theta deflection
cc = convertContact(C);
TH_d = nan(size(TH));
TH_d_alt =nan(size(TH));
for ii = 1:size(cc,1)
    
    idx = cc(ii,1);
    
    while isnan(TH(idx))
        idx = idx+1;
        if idx>=cc(ii,2) | idx >= length(TH)
            break
        end
    end
    
    x1 = xs{idx}(1);
    y1 = ys{idx}(1);
    % calculate angle between Contact point and basepoint at first point of
    % contact.
    TH_tc = atan2(CP(idx,2)-y1,CP(idx,1)-x1)*180/pi;
    for jj = 0:(cc(ii,2)-cc(ii,1))
        if length(xs{idx+jj})<1
            TH_d(idx+jj) = nan;
            continue
        end
        x1 = xs{idx+jj}(1);
        y1 = ys{idx+jj}(1);
        TH_d(idx+jj) = atan2(CP(idx+jj,2)-y1,CP(idx+jj,1)-x1)*180/pi - TH_tc;
    end
    
    
end

% wrap
TH_d(TH_d>nanmean(TH_d)+180) = TH_d(TH_d>nanmean(TH_d)+180)-360;
TH_d(TH_d<nanmean(TH_d)-180) = TH_d(TH_d<nanmean(TH_d)-180)+360;

E2D.TH_d = TH_d;
E2D_medfilt.TH_d = medfilt1(TH_d);


%% delta theta cp
TH_dcp = nan(size(TH));
for ii = 1:size(cc,1)
    idx = cc(ii,1);
    
    
    
    while isnan(TH(idx))
        idx = idx+1;
        if idx>cc(ii,2) | idx > length(TH)
            break
        end
    end
    
    
    TH_dcp(idx:cc(ii,2)) = TH(idx:cc(ii,2)) - TH(idx,1);
end

% wrap
TH_dcp(TH_dcp>nanmean(TH_dcp)+180) = TH_dcp(TH_dcp>nanmean(TH_dcp)+180)-360;
TH_dcp(TH_dcp<nanmean(TH_dcp)-180) = TH_dcp(TH_dcp<nanmean(TH_dcp)-180)+360;
% output
E2D.TH_dcp = TH_dcp;
E2D_medfilt.TH_dcp = medfilt1(TH_dcp);
%% Interp over nangaps

fNames = fieldnames(E2D);
for ii = 1:length(fNames)
    if ismember(fNames{ii},{'xs','ys','C','TH'})
        continue
    end
    if size(E2D.(fNames{ii}),2)>1
        tempUp = [];
        tempUpFilt = [];
        temp = E2D.(fNames{ii});
        tempFilt = E2D_medfilt.(fNames{ii});
        
        for jj = 1:size(temp,2)
            tempUp(:,jj) = InterpolateOverNans(temp(:,jj),20);
            tempUpFilt(:,jj) = InterpolateOverNans(tempFilt(:,jj),20);
        end
        E2D.(fNames{ii}) = tempUp;
        E2D_medfilt.(fNames{ii}) = tempUpFilt;
        continue
    end
    
    E2D_medfilt.(fNames{ii}) = InterpolateOverNans(E2D_medfilt.(fNames{ii}),20);
    E2D.(fNames{ii}) = InterpolateOverNans(E2D.(fNames{ii}),20);
end
%% Delete Outliers

fNames = fieldnames(E2D);
for ii = 1:length(fNames)
    if ismember(fNames{ii},{'xs','ys','BP','CP','TH','C'})
        E2D_upsamp.(fNames{ii}) = E2D.(fNames{ii});
        continue
    end
    
    E2D_medfilt.(fNames{ii}) = deleteoutliers(E2D_medfilt.(fNames{ii}),.001,1);
    
    E2D.(fNames{ii}) = deleteoutliers(E2D.(fNames{ii}),.0001,1);
end



%% Perform Upsampling
fNames = fieldnames(E2D);
for ii = 1:length(fNames)
    if ismember(fNames{ii},{'xs','ys','BP'})
        E2D_upsamp.(fNames{ii}) = E2D.(fNames{ii});
        continue
    end
    
    E2D_medfilt_upsamp.(fNames{ii}) = upsampForNeural(E2D_medfilt.(fNames{ii}),neural_word(:,1),frameCapTimes_ms,1,length(E2D.FX),length(E2D.FX));
    E2D_upsamp.(fNames{ii}) = upsampForNeural(E2D.(fNames{ii}),neural_word(:,1),frameCapTimes_ms,1,length(E2D.FX),length(E2D.FX));
end

%% Get R
R = nan(size(E2D_medfilt.C));
R_noFilt = R;
for ii = 1:length(E2D_medfilt.xs)
    if isempty(E2D_medfilt.xs{ii}) | isnan(E2D_medfilt.CP(ii,1))
        R(ii) = NaN;
        continue
    end
    R(ii) = sqrt((E2D_medfilt.xs{ii}(1)-E2D_medfilt.CP(ii,1)).^2 + (E2D_medfilt.ys{ii}(1)-E2D_medfilt.CP(ii,2)).^2);
    R_noFilt(ii) = sqrt((E2D.xs{ii}(1)-E2D.CP(ii,1)).^2 + (E2D.ys{ii}(1)-E2D.CP(ii,2)).^2);
end
%%
R = InterpolateOverNans(R,20);
R_noFilt = InterpolateOverNans(R_noFilt,20);

R_filt_upsamp = upsampForNeural(R,neural_word(:,1),frameCapTimes_ms,1,length(R),length(R));
R_upsamp = upsampForNeural(R_noFilt,neural_word(:,1),frameCapTimes_ms,1,length(R),length(R));

%% Upsample the RCCR logical
RCCR_upsamp = upsampForNeural(RCCR,neural_word(:,1),frameCapTimes_ms,1,length(E2D.C),length(E2D.C));
RCCR_upsamp(isnan(RCCR_upsamp)) = 0;
RCCR_upsamp = logical(RCCR_upsamp);
%% reformat into matrices with just the RCCR times
Mech.FX = E2D_upsamp.FX(RCCR_upsamp);
Mech.FY = E2D_upsamp.FY(RCCR_upsamp);
Mech.M = E2D_upsamp.M(RCCR_upsamp);

Mech.filtFX = E2D_medfilt_upsamp.FX(RCCR_upsamp);
Mech.filtFY = E2D_medfilt_upsamp.FY(RCCR_upsamp);
Mech.filtM = E2D_medfilt_upsamp.M(RCCR_upsamp);

Mech.filtAll = [Mech.filtFX Mech.filtFY Mech.filtM];

Geo.R = R_upsamp(RCCR_upsamp);
Geo.TH_cp = E2D_upsamp.TH_cp(RCCR_upsamp);
Geo.TH_d = E2D_upsamp.TH_d(RCCR_upsamp);
Geo.TH_dcp = E2D_upsamp.TH_dcp(RCCR_upsamp);


Geo.filtR = R_filt_upsamp(RCCR_upsamp);
Geo.filtTH_cp = E2D_medfilt_upsamp.TH_cp(RCCR_upsamp);
Geo.filtTH_d = E2D_medfilt_upsamp.TH_d(RCCR_upsamp);
Geo.filtTH_dcp = E2D_medfilt_upsamp.TH_dcp(RCCR_upsamp);


Geo.filtAll = [Geo.filtR Geo.filtTH_d];
C = E2D_upsamp.C(RCCR_upsamp);




%% Save output

overwrite = 1;
if exist([fOut '_cell_1*.mat'],'file')
    overwrite = input('File found, do you want to overwrite?(1/0)')
end

% [prox,med,dis] = getRadialDistanceGroup(Geo)

if overwrite
    
    for ii = 1:length(useCells)
        spikevec = neural_word(RCCR_upsamp,ii);
        neuralShapes = shapes(ii).shapes;
        save([fOut '_cell_' num2str(ii) '_toGLM.mat'],'Mech','Geo','spikevec','C')
    end
end

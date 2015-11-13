% alginNeuralAnd2D
clear
[fName,pName] = uigetfile('*varConcat.mat','Load in the concatenated variables from E2D');
load([pName fName]);

oldir = pwd;
cd C:\Users\nbush257\Documents\hartmann_lab\data\VG2D\neural\;
[fName,pName] = uigetfile('*.mat',['Load in the sorted neural data for ' fName]);
load([pName fName]);
cd(oldir);
fOut = regexp(fName,'_t\d{2}');fOut = fName(1:fOut+3);
%% Neural data to ms
frameCapTimes(frameCapTimes>time(length(time)))=NaN;
timeSubSamp = floor(time*1000);
timeMs = 1:timeSubSamp(end);
frameCapTimes_ms = frameCapTimes*1000;
figure
for ii = 1:length(shapes)
    subplot(1,length(shapes),ii)
    plot(shapes(ii).shapes','k')
    title(num2str(ii))
    axy(-150,150)
end

useCells = input('Which cells do you want to use?');
neural_clip = [];
neural_ms = [];
neural_word = zeros(length(timeMs),length(useCells));
for ii = 1:length(useCells)
    neural_clip = times(useCells(ii)).times; % pick your cell
    neural_ms = timeSubSamp(neural_clip);
    neural_word(neural_ms,ii) = 1;
end


%% getTH
xs = E2D.xs;
ys = E2D.ys;
C = E2D.C;
CP = E2D.CP;
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
if isrow(TH);TH = TH';end
E2D.TH = TH;
E2D_medfilt.TH = medfilt1(TH);

%% Interp over nangaps 

fNames = fieldnames(E2D);
for ii = 1:length(fNames)
    if ismember(fNames{ii},{'xs','ys','C'})
        E2D_upsamp.(fNames{ii}) = E2D.(fNames{ii});
        continue
    end
    if size(E2D.(fNames{ii}),2)>1
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


%% Perform Upsampling
fNames = fieldnames(E2D);
for ii = 1:length(fNames)
    if ismember(fNames{ii},{'xs','ys'})
        E2D_upsamp.(fNames{ii}) = E2D.(fNames{ii});
        continue
    end
    
    E2D_medfilt_upsamp.(fNames{ii}) = upsampForNeural(E2D_medfilt.(fNames{ii}),neural_word(:,1),frameCapTimes_ms,1,length(frameCapTimes),length(frameCapTimes));
    E2D_upsamp.(fNames{ii}) = upsampForNeural(E2D.(fNames{ii}),neural_word(:,1),frameCapTimes_ms,1,length(frameCapTimes),length(frameCapTimes));
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
R = InterpolateOverNans(R,20);
R_noFilt = InterpolateOverNans(R_noFilt,20);

R_filt_upsamp = upsampForNeural(R,neural_word(:,1),frameCapTimes_ms,1,length(frameCapTimes),length(frameCapTimes));
R_upsamp = upsampForNeural(R_noFilt,neural_word(:,1),frameCapTimes_ms,1,length(frameCapTimes),length(frameCapTimes));

%% Upsample the RCCR logical
RCCR_upsamp = upsampForNeural(RCCR,neural_word(:,1),frameCapTimes_ms,1,length(frameCapTimes),length(frameCapTimes));
RCCR_upsamp(isnan(RCCR_upsamp)) = 0;
RCCR_upsamp = logical(RCCR_upsamp);
%% reformat into matrices with just the RCCR times
Mech.FX = E2D_upsamp.FX(RCCR_upsamp);
Mech.FY = E2D_upsamp.FY(RCCR_upsamp);
Mech.M = E2D_upsamp.M(RCCR_upsamp);

Mech.filtFX = E2D_medfilt_upsamp.FX(RCCR_upsamp);
Mech.filtFY = E2D_medfilt_upsamp.FY(RCCR_upsamp);
Mech.filtM = E2D_medfilt_upsamp.M(RCCR_upsamp);

Geo.R = R_upsamp(RCCR_upsamp);
Geo.TH = E2D_upsamp.TH(RCCR_upsamp);

Geo.filtR = R_filt_upsamp(RCCR_upsamp);
Geo.filtTH = E2D_medfilt_upsamp.TH(RCCR_upsamp);



%% Save output

overwrite = 1;
if exist([fOut '_cell1 .mat'],'file')
    overwrite = input('File found, do you want to overwrite?(1/0)')
end

[prox,med,dis] = getRadialDistanceGroup(Geo);

if overwrite
    
    for ii = 1:length(useCells)
        spikevec = neural_word(:,ii);
        neuralShapes = shapes(ii).shapes;
        save([fOut '_cell' num2str(ii) '.mat'],'Mech','Geo','spikevec')
    end
end

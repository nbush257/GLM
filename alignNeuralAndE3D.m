function alignNeuralAndE3D
times = [];
[fName,pName] = uigetfile('*.mat','Load in the E3D data');
load([pName fName])
if any(strfind(fName,'E3D'))
    fOut = [pName fName(8:strfind(fName,'iter')-1) 'toGLM'];
else
    warning('Input filename does not contain ''E3D'', make sure this is the right file')
end
cd(pName)
[fName,pName] = uigetfile('*.mat',['Load in the sorted neural data for ' fName]);
load([pName fName])
%% Get timings
frameCapTimes(frameCapTimes>time(length(time)))=NaN;
timeSubSamp = floor(time*1000);
timeMs = 1:timeSubSamp(end);
frameCapTimes_ms = frameCapTimes*1000;
figure
for ii = 1:length(shapes)
    subplot(1,length(shapes),ii)
    plot(shapes(ii).shapes','k')
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
%% Median or LULU filtering for outlier removal
M = medfilt1(M,3);
F  =medfilt1(F,3);
PHIE_filt = medfilt1(PHIE_filt);
TH_filt = medfilt1(TH_filt);
Rcp = medfilt1(Rcp);

%% upsample
for ii = 1:3
    M_upsamp(:,ii) = upsampForNeural(M(:,ii),neural_word(:,1),frameCapTimes_ms,1,length(frameCapTimes),length(frameCapTimes));
    F_upsamp(:,ii) = upsampForNeural(F(:,ii),neural_word(:,1),frameCapTimes_ms,1,length(frameCapTimes),length(frameCapTimes));
end

PHIE_upsamp = upsampForNeural(PHIE_filt,neural_word(:,1),frameCapTimes_ms,1,length(frameCapTimes),length(frameCapTimes));
Rupsamp= upsampForNeural(Rcp,neural_word(:,1),frameCapTimes_ms,1,length(frameCapTimes),length(frameCapTimes));
TH_upsamp = upsampForNeural(TH_filt,neural_word(:,1),frameCapTimes_ms,1,length(frameCapTimes),length(frameCapTimes));
ZETA_upsamp = upsampForNeural(ZETA_filt,neural_word(:,1),frameCapTimes_ms,1,length(frameCapTimes),length(frameCapTimes));
%% Data out
Mech.Mx = M_upsamp(:,1);
Mech.My = M_upsamp(:,2);
Mech.Mz = M_upsamp(:,3);

Mech.Fx = F_upsamp(:,1);
Mech.Fy = F_upsamp(:,2);
Mech.Fz = F_upsamp(:,3);

Mech.all = [Mech.Mx Mech.My Mech.Mz Mech.Fx Mech.Fy Mech.Fz];
Mech.order = ['Mx','My','Mz','Fx','Fy','Fz'];

Geo.R = Rupsamp;
Geo.TH = TH_upsamp;
Geo.Zeta = ZETA_upsamp;
Geo.Phi = PHIE_upsamp;

Geo.all = [Geo.R Geo.TH Geo.Zeta Geo.Phi];
Geo.order = ['R','TH','Zeta','Phi'];
C = zeros(size(Geo.R));
C(~isnan(Geo.R))=1;
%% Save
overwrite = 1;
if exist([fOut '_cell1 .mat'],'file')
    overwrite = input('File found, do you want to overwrite?(1/0)')
end

[prox,med,dis] = getRadialDistanceGroup(Geo);

if overwrite
    
    for ii = 1:length(useCells)
        neuralOut = neural_word(:,ii);
        neuralShapes = shapes(ii).shapes;
        save([fOut '_cell' num2str(ii) '.mat'],'Mech','Geo','neuralOut','neuralShapes','prox','med','dis')
    end
end




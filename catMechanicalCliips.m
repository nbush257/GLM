clear; ca
whisker = 'C1'
ratNum = '08'
trial = 't01'
startFrame = 1;
cell = 1;
sr = 40000;
RCCR = [5500 17958];


%%
dProc = dir(['*PROC*' whisker '*' trial '*']);
allFX = [];
allFY = [];
allM = [];
allBP = [];
allC = [];
allTH = [];
allxs = [];
allys = [];
allCP = [];

dE2D = dir(['*E2D*' whisker '*' trial '*']);
for ii = 1:length(dE2D)
load(dE2D(ii).name)
load(dProc(ii).name)
% if ii ==2
%     FX = nan(20000,1);
%     FY = nan(20000,1);
%     M = nan(20000,1);
%     BP = nan(20000,2);
%     C = zeros(20000,1);
%     TH = nan(20000,1);
%     xs{20000} = [];
%     ys{20000} = [];
%     CP = nan(20000,2);
% end
%     
    
    
allFX = [allFX;FX];
allFY = [allFY;FY];
allM = [allM;M];
allBP = [allBP;BP];
allC = [allC;C];
allTH = [allTH;TH];
allxs = [allxs xs];
allys = [allys ys];
allCP = [allCP;CP];
end

dNeural = dir(['C:\Users\nbush257\Documents\hartmann_lab\data\VG2D\neural\*' ratNum '*' whisker '*' trial '*sorted*']);
load(['C:\Users\nbush257\Documents\hartmann_lab\data\VG2D\neural\' dNeural.name])
disp('===============================')
disp('')

disp(['You are on cell ' num2str(cell) ' of ' num2str(length(times))])

endFrame = length(allC)+startFrame-1;

FX = allFX;
FY = allFY;
M = allM;
BP = allBP;
C = allC;
TH = allTH;
xs = allxs;
ys = allys;
CP = allCP;
TAG = dNeural.name(1:end-11);
%%
% keyboard
glmPrep_1khz;
newRCCR = [5500 17958;40001+13375 40001+19999;60001 60001+2000;80001+12000 80001+19999];
temp = zeros(size(mechMat_25,1),1);
for ii=1:4
tStart= frameCapSamps(newRCCR(ii,1)+startFrame);
tEnd = frameCapSamps(newRCCR(ii,2)+startFrame);

RCCR_start_ms = timeSubSamp(tStart);
RCCR_end_ms = timeSubSamp(tEnd);
temp(RCCR_start_ms:RCCR_end_ms) = 1;
end
mech_25_justRCCR = mechMat_25(temp==1,:);
mech_85_justRCCR = mechMat_25(temp==1,:);
mech_noFilt_justRCCR = mechMat_25(temp==1,:);
geo_25_justRCCR = mechMat_25(temp==1,:);
geo_85_justRCCR = mechMat_25(temp==1,:);
geo_noFilt_justRCCR = mechMat_25(temp==1,:);

save([TAG '_cell_' num2str(cell) '_merged.mat'],'geo_25_justRCCR','geo_85_justRCCR','geo_noFilt_justRCCR','mech_25_justRCCR','mech_85_justRCCR','mech_noFilt_justRCCR','C_upsamp','shapes','times')


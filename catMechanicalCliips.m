clear; ca
whisker = 'D4'
trial = 't01'
startFrame = 1;
cell = 1;
sr = 40000;
RCCR = [1 60000 + 11800];


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

dNeural = dir(['D:\data\analyzed\2015_15\*' whisker '*' trial '*sorted*']);
load(['D:\data\analyzed\2015_15\' dNeural.name])
disp('===============================')
disp('')

disp(['You are on cell ' num2str(cell) ' of ' num2str(length(times))])

endFrame = length(allC);

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
glmPrep_1khz;
save([TAG '_cell_' num2str(cell) '_merged.mat'],'geo_25_justRCCR','geo_85_justRCCR','geo_noFilt_justRCCR','mech_25_justRCCR','mech_85_justRCCR','mech_noFilt_justRCCR','C_upsamp','shapes','times')

